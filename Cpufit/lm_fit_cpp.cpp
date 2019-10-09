#include "cpufit.h"
#include "../Gpufit/constants.h"
#include "lm_fit.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>

// TODO if std::size_t and int are not the same, we will get lots of C26451 warnings here related to it, they can be ignored or
// int should be converted to size_t but be careful, there is at least one for loop that checks for >=0 which only works with int that way
// MS C compiler 16.1 (2019) shows the behavior for example

LMFitCPP::LMFitCPP(
    REAL const tolerance,
    std::size_t const fit_index,
    REAL const * data,
    REAL const * weight,
    Info const & info,
    REAL const * initial_parameters,
    int const * parameters_to_fit,
    char * user_info,
    REAL * output_parameters,
    int * output_state,
    REAL * output_chi_square,
    int * output_n_iterations
    ) :
    fit_index_(fit_index),
    data_(data),
    weight_(weight),
    initial_parameters_(initial_parameters),
    tolerance_(tolerance),
    converged_(false),
    info_(info),
    parameters_to_fit_(parameters_to_fit),
    curve_(info.n_points_),
    derivatives_(info.n_points_*info.n_parameters_),
    hessian_(info.n_parameters_to_fit_*info.n_parameters_to_fit_),
    modified_hessian_(info.n_parameters_to_fit_*info.n_parameters_to_fit_),
    decomposed_hessian_(info.n_parameters_to_fit_*info.n_parameters_to_fit_),
    pivot_array_(info_.n_parameters_to_fit_),
    gradient_(info.n_parameters_to_fit_),
    delta_(info.n_parameters_to_fit_),
    scaling_vector_(info.n_parameters_to_fit_),
    prev_chi_square_(0),
    lambda_(0.001f),
    prev_parameters_(info.n_parameters_to_fit_),
    user_info_(user_info),
    parameters_(output_parameters),
    state_(output_state),
    chi_square_(output_chi_square),
    n_iterations_(output_n_iterations)
{}

template<class T>
int decompose_LUP(std::vector<T> & matrix, int const N, double const Tol, std::vector<int> & permutation_vector) {

    for (int i = 0; i < N; i++)
        permutation_vector[i] = i;

    for (int i = 0; i < N; i++)
    {
        T max_value = 0;
        int max_index = i;

        for (int k = i; k < N; k++)
        {
            T absolute_value = std::abs(matrix[k * N + i]);
            if (absolute_value > max_value)
            {
                max_value = absolute_value;
                max_index = k;
            }
        }

        if (max_value < Tol)
            return 0; //failure, matrix is degenerate

        if (max_index != i)
        {
            //pivoting permutation vector
            std::swap(permutation_vector[i], permutation_vector[max_index]);

            //pivoting rows of matrix
            for (int j = 0; j < N; j++)
                std::swap(matrix[i * N + j], matrix[max_index * N + j]);
        }

        for (int j = i + 1; j < N; j++)
        {
            matrix[j * N + i] /= matrix[i * N + i];

            for (int k = i + 1; k < N; k++)
                matrix[j * N + k] -= matrix[j * N + i] * matrix[i * N + k];
        }
    }

    return 1;  //decomposition done 
}

template<class T>
void solve_LUP(
    std::vector<T> const & matrix,
    std::vector<int> const & permutation_vector,
    std::vector<T> const & vector,
    int const N,
    std::vector<T> & solution)
{
    for (int i = 0; i < N; i++)
    {
        solution[i] = vector[permutation_vector[i]];

        for (int k = 0; k < i; k++)
        {
            solution[i] -= matrix[i * N + k] * solution[k];
        }
    }

    for (int i = N - 1; i >= 0; i--)
    {
        for (int k = i + 1; k < N; k++)
        {
            solution[i] -= matrix[i * N + k] * solution[k];
        }

        solution[i] = solution[i] / matrix[i * N + i];
    }
}

void LMFitCPP::decompose_hessian_LUP(std::vector<REAL> const & hessian)
{
    decomposed_hessian_ = hessian;

    int const N = int(gradient_.size());

    int const singular = decompose_LUP(decomposed_hessian_, info_.n_parameters_to_fit_, 0.0, pivot_array_);
    if (singular == 0)
        *state_ = FitState::SINGULAR_HESSIAN;
}

void LMFitCPP::calc_derivatives_gauss2d(
    std::vector<REAL> & derivatives)
{
    std::size_t const  fit_size_x = std::size_t(std::sqrt(info_.n_points_));

    for (std::size_t y = 0; y < fit_size_x; y++)
        for (std::size_t x = 0; x < fit_size_x; x++)
        {
            REAL const argx = (x - parameters_[1]) * (x - parameters_[1]) / (2 * parameters_[3] * parameters_[3]);
            REAL const argy = (y - parameters_[2]) * (y - parameters_[2]) / (2 * parameters_[3] * parameters_[3]);
            REAL const ex = exp(-(argx + argy));

            derivatives[0 * info_.n_points_ + y*fit_size_x + x]
                = ex;
            derivatives[1 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (x - parameters_[1])*ex) / (parameters_[3] * parameters_[3]);
            derivatives[2 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (y - parameters_[2])*ex) / (parameters_[3] * parameters_[3]);
            derivatives[3 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0]
                * ((x - parameters_[1])*(x - parameters_[1])
                + (y - parameters_[2])*(y - parameters_[2]))*ex)
                / (parameters_[3] * parameters_[3] * parameters_[3]);
            derivatives[4 * info_.n_points_ + y*fit_size_x + x]
                = 1;
        }
}

void LMFitCPP::calc_derivatives_gauss2delliptic(
    std::vector<REAL> & derivatives)
{
    std::size_t const  fit_size_x = std::size_t(std::sqrt(info_.n_points_));

    for (std::size_t y = 0; y < fit_size_x; y++)
        for (std::size_t x = 0; x < fit_size_x; x++)
        {
            REAL const argx = (x - parameters_[1]) * (x - parameters_[1]) / (2 * parameters_[3] * parameters_[3]);
            REAL const argy = (y - parameters_[2]) * (y - parameters_[2]) / (2 * parameters_[4] * parameters_[4]);
            REAL const ex = exp(-(argx +argy));

            derivatives[0 * info_.n_points_ + y*fit_size_x + x]
                = ex;
            derivatives[1 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (x - parameters_[1])*ex) / (parameters_[3] * parameters_[3]);
            derivatives[2 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (y - parameters_[2])*ex) / (parameters_[4] * parameters_[4]);
            derivatives[3 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (x - parameters_[1])*(x - parameters_[1])*ex) / (parameters_[3] * parameters_[3] * parameters_[3]);
            derivatives[4 * info_.n_points_ + y*fit_size_x + x]
                = (parameters_[0] * (y - parameters_[2])*(y - parameters_[2])*ex) / (parameters_[4] * parameters_[4] * parameters_[4]);
            derivatives[5 * info_.n_points_ + y*fit_size_x + x]
                = 1;
        }
}

void LMFitCPP::calc_derivatives_gauss2drotated(
    std::vector<REAL> & derivatives)
{
    std::size_t const  fit_size_x = std::size_t(std::sqrt(info_.n_points_));

    REAL const amplitude = parameters_[0];
    REAL const x0 = parameters_[1];
    REAL const y0 = parameters_[2];
    REAL const sig_x = parameters_[3];
    REAL const sig_y = parameters_[4];
    REAL const background = parameters_[5];
    REAL const rot_sin = sin(parameters_[6]);
    REAL const rot_cos = cos(parameters_[6]);

    for (std::size_t y = 0; y < fit_size_x; y++)
        for (std::size_t x = 0; x < fit_size_x; x++)
        {
            REAL const arga = ((x - x0) * rot_cos) - ((y - y0) * rot_sin);
            REAL const argb = ((x - x0) * rot_sin) + ((y - y0) * rot_cos);
            REAL const ex = exp((-0.5f) * (((arga / sig_x) * (arga / sig_x)) + ((argb / sig_y) * (argb / sig_y))));

            derivatives[0 * info_.n_points_ + y*fit_size_x + x]
                = ex;
            derivatives[1 * info_.n_points_ + y*fit_size_x + x]
                = ex * (amplitude * rot_cos * arga / (sig_x*sig_x) + amplitude * rot_sin *argb / (sig_y*sig_y));
            derivatives[2 * info_.n_points_ + y*fit_size_x + x]
                = ex * (-amplitude * rot_sin * arga / (sig_x*sig_x) + amplitude * rot_cos *argb / (sig_y*sig_y));
            derivatives[3 * info_.n_points_ + y*fit_size_x + x]
                = ex * amplitude * arga * arga / (sig_x*sig_x*sig_x);
            derivatives[4 * info_.n_points_ + y*fit_size_x + x]
                = ex * amplitude * argb * argb / (sig_y*sig_y*sig_y);
            derivatives[5 * info_.n_points_ + y*fit_size_x + x]
                = 1.f;
            derivatives[6 * info_.n_points_ + y*fit_size_x + x]
                = ex * amplitude * arga * argb * (1.f / (sig_x*sig_x) - 1.f / (sig_y*sig_y));
        }
}

void LMFitCPP::calc_derivatives_gauss1d(
    std::vector<REAL> & derivatives)
{
    REAL * user_info_float = (REAL*)user_info_;
    REAL x = 0.;

    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        if (!user_info_float)
        {
            x = REAL(point_index);
        }
        else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
        {
            x = user_info_float[point_index];
        }
        else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_)
        {
            std::size_t const fit_begin = fit_index_ * info_.n_points_;
            x = user_info_float[fit_begin + point_index];
        }

        REAL argx = ((x - parameters_[1])*(x - parameters_[1])) / (2 * parameters_[2] * parameters_[2]);
        REAL ex = exp(-argx);

        derivatives[0 * info_.n_points_ + point_index] = ex;
        derivatives[1 * info_.n_points_ + point_index] = (parameters_[0] * (x - parameters_[1])*ex) / (parameters_[2] * parameters_[2]);
        derivatives[2 * info_.n_points_ + point_index] = (parameters_[0] * (x - parameters_[1])*(x - parameters_[1])*ex) / (parameters_[2] * parameters_[2] * parameters_[2]);
        derivatives[3 * info_.n_points_ + point_index] = 1;
    }
}

void LMFitCPP::calc_derivatives_exp1d(std::vector<REAL>& derivatives){
  REAL * user_info_float = (REAL*)user_info_;
  REAL x = 0.;
  for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++){
    if (!user_info_float)
      x = REAL(point_index);
    else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
      x = user_info_float[point_index];
    else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_){
      std::size_t const fit_begin = fit_index_ * info_.n_points_;
      x = user_info_float[fit_begin + point_index];
    }
    REAL ev = exp(-x*parameters_[2]);
    derivatives[0*info_.n_points_+point_index] = 1;
    derivatives[1*info_.n_points_+point_index] = ev;
    derivatives[2*info_.n_points_+point_index] = -x*ev;
  }
}

void LMFitCPP::calc_derivatives_erfc_xgauss(std::vector<REAL>& derivatives){
  REAL * user_info_float = (REAL*)user_info_;
  REAL x = 0.;
  REAL k2 = parameters_[1] / parameters_[2];
  REAL k3 = k2 * k2;
  for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++){
    if (!user_info_float)
      x = REAL(point_index);
    else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
      x = user_info_float[point_index];
    else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_){
      std::size_t const fit_begin = fit_index_ * info_.n_points_;
      x = user_info_float[fit_begin + point_index];
    }
    REAL k0 = (x-parameters_[0]) / parameters_[1];
    REAL k1 = (x-parameters_[0]) / parameters_[2];
    REAL v0 = 0.5f / abs(parameters_[2]);
    v0 *= expf(k1-0.5f*k3);
    v0 *= erfcf(parameters_[2]/abs(parameters_[2])/sqrtf(2.f)*(k0+k2));
    REAL v1 = v0 / parameters_[2];
    REAL ev=expf(-k3-0.5f*k0*k0)/parameters_[1]/parameters_[2]/sqrtf(2.f*3.14159f);
    derivatives[0*info_.n_points_+point_index] = parameters_[3]*(-v1*k2+(k0-k2)*ev);
    derivatives[1*info_.n_points_+point_index] = parameters_[3]*(-v1+ev);
    derivatives[2*info_.n_points_+point_index] = parameters_[3]*(v1*(k3-1)+k3*ev);
    derivatives[3*info_.n_points_+point_index] = v0;
    derivatives[4*info_.n_points_+point_index] = 1.0f;
  }
}

void LMFitCPP::calc_derivatives_cauchy2delliptic(
    std::vector<REAL> & derivatives)
{
    std::size_t const  fit_size_x = std::size_t(std::sqrt(info_.n_points_));

    for (std::size_t y = 0; y < fit_size_x; y++)
        for (std::size_t x = 0; x < fit_size_x; x++)
        {
            REAL const argx =
                ((parameters_[1] - x) / parameters_[3])
                *((parameters_[1] - x) / parameters_[3]) + 1.f;
            REAL const argy =
                ((parameters_[2] - y) / parameters_[4])
                *((parameters_[2] - y) / parameters_[4]) + 1.f;

            derivatives[0 * info_.n_points_ + y*fit_size_x + x]
                = 1.f / (argx*argy);
            derivatives[1 * info_.n_points_ + y*fit_size_x + x] =
                -2.f * parameters_[0] * (parameters_[1] - x)
                / (parameters_[3] * parameters_[3] * argx*argx*argy);
            derivatives[2 * info_.n_points_ + y*fit_size_x + x] =
                -2.f * parameters_[0] * (parameters_[2] - y)
                / (parameters_[4] * parameters_[4] * argy*argy*argx);
            derivatives[3 * info_.n_points_ + y*fit_size_x + x] =
                2.f * parameters_[0] * (parameters_[1] - x) * (parameters_[1] - x)
                / (parameters_[3] * parameters_[3] * parameters_[3] * argx*argx*argy);
            derivatives[4 * info_.n_points_ + y*fit_size_x + x] =
                2.f * parameters_[0] * (parameters_[2] - y) * (parameters_[2] - y)
                / (parameters_[4] * parameters_[4] * parameters_[4] * argy*argy*argx);
            derivatives[5 * info_.n_points_ + y*fit_size_x + x]
                = 1.f;
        }
}

void LMFitCPP::calc_derivatives_linear1d(
    std::vector<REAL> & derivatives)
{
    REAL * user_info_float = (REAL*)user_info_;
    REAL x = 0.;

    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        if (!user_info_float)
        {
            x = REAL(point_index);
        }
        else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
        {
            x = user_info_float[point_index];
        }
        else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_)
        {
            std::size_t const fit_begin = fit_index_ * info_.n_points_;
            x = user_info_float[fit_begin + point_index];
        }

        derivatives[0 * info_.n_points_ + point_index] = 1.;
        derivatives[1 * info_.n_points_ + point_index] = x;
    }
}

void LMFitCPP::calc_derivatives_fletcher_powell_helix(
    std::vector<REAL> & derivatives)
{
    REAL const pi = 3.14159f;

    REAL const * p = parameters_;

    REAL const arg = p[0] * p[0] + p[1] * p[1];

    // derivatives with respect to p[0]
    derivatives[0 * info_.n_points_ + 0] = 100.f * 1.f / (2.f*pi) * p[1] / arg;
    derivatives[0 * info_.n_points_ + 1] = 10.f * p[0] / std::sqrt(arg);
    derivatives[0 * info_.n_points_ + 2] = 0.f;

    // derivatives with respect to p[1]
    derivatives[1 * info_.n_points_ + 0] = -100.f * 1.f / (2.f*pi) * p[0] / (arg);
    derivatives[1 * info_.n_points_ + 1] = 10.f * p[1] / std::sqrt(arg);
    derivatives[1 * info_.n_points_ + 2] = 0.f;

    // derivatives with respect to p[2]
    derivatives[2 * info_.n_points_ + 0] = 10.f;
    derivatives[2 * info_.n_points_ + 1] = 0.f;
    derivatives[2 * info_.n_points_ + 2] = 1.f;
}

void LMFitCPP::calc_derivatives_brown_dennis(
    std::vector<REAL> & derivatives)
{
    REAL const * p = parameters_;

    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        REAL const t = static_cast<REAL>(point_index) / 5.f;

        REAL const arg1 = p[0] + p[1] * t - std::exp(t);
        REAL const arg2 = p[2] + p[3] * std::sin(t) - std::cos(t);

        derivatives[0 * info_.n_points_ + point_index] = 2.f * arg1;
        derivatives[1 * info_.n_points_ + point_index] = 2.f * t * arg1;
        derivatives[2 * info_.n_points_ + point_index] = 2.f * arg2;
        derivatives[3 * info_.n_points_ + point_index] = 2.f * std::sin(t) * arg2;
    }
}

void LMFitCPP::calc_values_cauchy2delliptic(std::vector<REAL>& cauchy)
{
    int const size_x = int(std::sqrt(REAL(info_.n_points_)));
    int const size_y = size_x;

    for (int iy = 0; iy < size_y; iy++)
    {
        for (int ix = 0; ix < size_x; ix++)
        {
            REAL const argx =
                ((parameters_[1] - ix) / parameters_[3])
                *((parameters_[1] - ix) / parameters_[3]) + 1.f;
            REAL const argy =
                ((parameters_[2] - iy) / parameters_[4])
                *((parameters_[2] - iy) / parameters_[4]) + 1.f;

            cauchy[iy*size_x + ix] = parameters_[0] / (argx * argy) + parameters_[5];
        }
    }
}

void LMFitCPP::calc_values_gauss2d(std::vector<REAL>& gaussian)
{
    int const size_x = int(std::sqrt(REAL(info_.n_points_)));
    int const size_y = size_x;

    for (int iy = 0; iy < size_y; iy++)
    {
        for (int ix = 0; ix < size_x; ix++)
        {
            REAL argx = (ix - parameters_[1]) * (ix - parameters_[1]) / (2 * parameters_[3] * parameters_[3]);
            REAL argy = (iy - parameters_[2]) * (iy - parameters_[2]) / (2 * parameters_[3] * parameters_[3]);
            REAL ex = exp(-(argx +argy));

            gaussian[iy*size_x + ix] = parameters_[0] * ex + parameters_[4];
        }
    }
}

void LMFitCPP::calc_values_gauss2delliptic(std::vector<REAL>& gaussian)
{
    int const size_x = int(std::sqrt(REAL(info_.n_points_)));
    int const size_y = size_x;
    for (int iy = 0; iy < size_y; iy++)
    {
        for (int ix = 0; ix < size_x; ix++)
        {
            REAL argx = (ix - parameters_[1]) * (ix - parameters_[1]) / (2 * parameters_[3] * parameters_[3]);
            REAL argy = (iy - parameters_[2]) * (iy - parameters_[2]) / (2 * parameters_[4] * parameters_[4]);
            REAL ex = exp(-(argx + argy));

            gaussian[iy*size_x + ix]
                = parameters_[0] * ex + parameters_[5];
        }
    }
}
    
void LMFitCPP::calc_values_gauss2drotated(std::vector<REAL>& gaussian)
{
    int const size_x = int(std::sqrt(REAL(info_.n_points_)));
    int const size_y = size_x;

    REAL amplitude = parameters_[0];
    REAL background = parameters_[5];
    REAL x0 = parameters_[1];
    REAL y0 = parameters_[2];
    REAL sig_x = parameters_[3];
    REAL sig_y = parameters_[4];
    REAL rot_sin = sin(parameters_[6]);
    REAL rot_cos = cos(parameters_[6]);

    for (int iy = 0; iy < size_y; iy++)
    {
        for (int ix = 0; ix < size_x; ix++)
        {
            int const pixel_index = iy*size_x + ix;

            REAL arga = ((ix - x0) * rot_cos) - ((iy - y0) * rot_sin);
            REAL argb = ((ix - x0) * rot_sin) + ((iy - y0) * rot_cos);

            REAL ex
                = exp((-0.5f) * (((arga / sig_x) * (arga / sig_x)) + ((argb / sig_y) * (argb / sig_y))));

            gaussian[pixel_index] = amplitude * ex + background;
        }
    }
}

void LMFitCPP::calc_values_gauss1d(std::vector<REAL>& gaussian)
{
    REAL * user_info_float = (REAL*)user_info_;
    REAL x = 0.f;
    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        if (!user_info_float)
        {
            x = REAL(point_index);
        }
        else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
        {
            x = user_info_float[point_index];
        }
        else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_)
        {
            std::size_t const fit_begin = fit_index_ * info_.n_points_;
            x = user_info_float[fit_begin + point_index];
        }

        REAL argx
            = ((x - parameters_[1])*(x - parameters_[1]))
            / (2.f * parameters_[2] * parameters_[2]);
        REAL ex = exp(-argx);
        gaussian[point_index] = parameters_[0] * ex + parameters_[3];
    }
}

void LMFitCPP::calc_values_exp1d(std::vector<REAL>& exponential)
{
  REAL * user_info_float = (REAL*)user_info_;
  REAL x = 0.f;
  for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++){
    if (!user_info_float)
      x = REAL(point_index);
    else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
      x = user_info_float[point_index];
    else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_){
      std::size_t const fit_begin = fit_index_ * info_.n_points_;
      x = user_info_float[fit_begin + point_index];
    }
    exponential[point_index] = parameters_[0]+parameters_[1]*exp(-x*parameters_[2]);
  }
}

void LMFitCPP::calc_values_erfc_xgauss(std::vector<REAL>& xgauss){
  REAL * user_info_float = (REAL*)user_info_;
  REAL x = 0.f;
  REAL k2 = parameters_[1] / parameters_[2];
  REAL k3 = k2 * k2;
  for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++){
    if (!user_info_float)
      x = REAL(point_index);
    else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
      x = user_info_float[point_index];
    else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_){
      std::size_t const fit_begin = fit_index_ * info_.n_points_;
      x = user_info_float[fit_begin + point_index];
    }
    REAL k0 = (x-parameters_[0]) / parameters_[1];
    REAL k1 = (x-parameters_[0]) / parameters_[2];    
    REAL v0 = 0.5f / abs(parameters_[2]);
    v0 *= expf(k1-0.5f*k3);
    v0 *= erfcf(parameters_[2]/abs(parameters_[2])/sqrtf(2.f)*(k0+k2));
    xgauss[point_index] = parameters_[3]*v0+parameters_[4];
  }
}
  
void LMFitCPP::calc_values_linear1d(std::vector<REAL>& line)
{
    REAL * user_info_float = (REAL*)user_info_;
    REAL x = 0.f;
    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        if (!user_info_float)
        {
            x = REAL(point_index);
        }
        else if (info_.user_info_size_ / sizeof(REAL) == info_.n_points_)
        {
            x = user_info_float[point_index];
        }
        else if (info_.user_info_size_ / sizeof(REAL) > info_.n_points_)
        {
            std::size_t const fit_begin = fit_index_ * info_.n_points_;
            x = user_info_float[fit_begin + point_index];
        }
        line[point_index] = parameters_[0] + parameters_[1] * x;
    }
}

void LMFitCPP::calc_values_fletcher_powell_helix(std::vector<REAL>& values)
{
    REAL const * p = parameters_;

    REAL const pi = 3.14159f;

    REAL theta = 0.f;

    if (0. < p[0])
        theta = .5f * atan(p[1] / p[0]) / pi;
    else if (p[0] < 0.)
        theta = .5f * atan(p[1] / p[0]) / pi + .5f;
    else if (0. < p[1])
        theta = .25f;
    else if (p[1] < 0.)
        theta = -.25f;
    else
        theta = 0.f;

    values[0] = 10.f * (p[2] - 10.f * theta);
    values[1] = 10.f * (std::sqrt(p[0] * p[0] + p[1] * p[1]) - 1.f);
    values[2] = p[2];
}

void LMFitCPP::calc_values_brown_dennis(std::vector<REAL>& values)
{
    REAL const * p = parameters_;

    for (std::size_t point_index = 0; point_index < info_.n_points_; point_index++)
    {
        REAL const t = static_cast<REAL>(point_index) / 5.f;

        REAL const arg1 = p[0] + p[1] * t - std::exp(t);
        REAL const arg2 = p[2] + p[3] * std::sin(t) - std::cos(t);

        values[point_index] = arg1*arg1 + arg2*arg2;
    }
}

void LMFitCPP::calc_curve_values(std::vector<REAL>& curve, std::vector<REAL>& derivatives)
{           
    if (info_.model_id_ == GAUSS_1D)
    {
        calc_values_gauss1d(curve);
        calc_derivatives_gauss1d(derivatives);
    }
    else if (info_.model_id_ == GAUSS_2D)
    {
        calc_values_gauss2d(curve);
        calc_derivatives_gauss2d(derivatives);
    }
    else if (info_.model_id_ == GAUSS_2D_ELLIPTIC)
    {
        calc_values_gauss2delliptic(curve);
        calc_derivatives_gauss2delliptic(derivatives);
    }
    else if (info_.model_id_ == GAUSS_2D_ROTATED)
    {
        calc_values_gauss2drotated(curve);
        calc_derivatives_gauss2drotated(derivatives);
    }
    else if (info_.model_id_ == CAUCHY_2D_ELLIPTIC)
    {
        calc_values_cauchy2delliptic(curve);
        calc_derivatives_cauchy2delliptic(derivatives);
    }
    else if (info_.model_id_ == LINEAR_1D)
    {
        calc_values_linear1d(curve);
        calc_derivatives_linear1d(derivatives);
    }
    else if (info_.model_id_ == EXP_1D)
    {
      calc_values_exp1d(curve);
      calc_derivatives_exp1d(derivatives);
    }
    else if (info_.model_id_ == ERFC_XGAUSS)
    {
      calc_values_erfc_xgauss(curve);
      calc_derivatives_erfc_xgauss(derivatives);
    }
    else if (info_.model_id_ == FLETCHER_POWELL_HELIX)
    {
        calc_values_fletcher_powell_helix(curve);
        calc_derivatives_fletcher_powell_helix(derivatives);
    }
    else if (info_.model_id_ == BROWN_DENNIS)
    {
        calc_values_brown_dennis(curve);
        calc_derivatives_brown_dennis(derivatives);
    }
}

void LMFitCPP::calculate_hessian(
    std::vector<REAL> const & derivatives,
    std::vector<REAL> const & curve)
{
    for (int jp = 0, jhessian = 0; jp < info_.n_parameters_; jp++)
    {
        if (parameters_to_fit_[jp])
        {
            for (int ip = 0, ihessian = 0; ip < jp + 1; ip++)
            {
                if (parameters_to_fit_[ip])
                {
                    std::size_t const ijhessian
                        = ihessian * info_.n_parameters_to_fit_ + jhessian;
                    std::size_t const jihessian
                        = jhessian * info_.n_parameters_to_fit_ + ihessian;
                    std::size_t const derivatives_index_i = ip*info_.n_points_;
                    std::size_t const derivatives_index_j = jp*info_.n_points_;
                    
                    double sum = 0.0;
                    for (std::size_t pixel_index = 0; pixel_index < info_.n_points_; pixel_index++)
                    {
                        if (info_.estimator_id_ == LSE)
                        {
                            if (!weight_)
                            {
                                sum
                                    += derivatives[derivatives_index_i + pixel_index]
                                    * derivatives[derivatives_index_j + pixel_index];
                            }
                            else
                            {
                                sum
                                    += derivatives[derivatives_index_i + pixel_index]
                                    * derivatives[derivatives_index_j + pixel_index]
                                    * weight_[pixel_index];
                            }
                        }
                        else if (info_.estimator_id_ == MLE)
                        {
                            sum
                                += data_[pixel_index] / (curve[pixel_index] * curve[pixel_index])
                                * derivatives[derivatives_index_i + pixel_index]
                                * derivatives[derivatives_index_j + pixel_index];
                        }
                    }
                    hessian_[ijhessian] = REAL(sum);
                    if (ijhessian != jihessian)
                    {
                        hessian_[jihessian]
                            = hessian_[ijhessian];
                    }
                    ihessian++;
                }
            }
            jhessian++;
        }
    }

}

void LMFitCPP::calc_gradient(
    std::vector<REAL> const & derivatives,
    std::vector<REAL> const & curve)
{

    for (int ip = 0, gradient_index = 0; ip < info_.n_parameters_; ip++)
    {
        if (parameters_to_fit_[ip])
        {
            std::size_t const derivatives_index = ip*info_.n_points_;
            double sum = 0.;
            for (std::size_t pixel_index = 0; pixel_index < info_.n_points_; pixel_index++)
            {
                REAL deviant = data_[pixel_index] - curve[pixel_index];

                if (info_.estimator_id_ == LSE)
                {
                    if (!weight_)
                    {
                        sum
                            += deviant * derivatives[derivatives_index + pixel_index];
                    }
                    else
                    {
                        sum
                            += deviant * derivatives[derivatives_index + pixel_index] * weight_[pixel_index];
                    }

                }
                else if (info_.estimator_id_ == MLE)
                {
                    sum
                        += -derivatives[derivatives_index + pixel_index] * (1 - data_[pixel_index] / curve[pixel_index]);
                }
            }
            gradient_[gradient_index] = REAL(sum);
            gradient_index++;
        }
    }

}

void LMFitCPP::calc_chi_square(
    std::vector<REAL> const & values)
{
    double sum = 0.0;
    for (size_t pixel_index = 0; pixel_index < values.size(); pixel_index++)
    {
        REAL deviant = values[pixel_index] - data_[pixel_index];
        if (info_.estimator_id_ == LSE)
        {
            if (!weight_)
            {
                sum += deviant * deviant;
            }
            else
            {
                sum += deviant * deviant * weight_[pixel_index];
            }
        }
        else if (info_.estimator_id_ == MLE)
        {
            if (values[pixel_index] <= 0.f)
            {
                *state_ = FitState::NEG_CURVATURE_MLE;
                return;
            }
            if (data_[pixel_index] != 0.f)
            {
                sum
                    += 2 * (deviant - data_[pixel_index] * std::log(values[pixel_index] / data_[pixel_index]));
            }
            else
            {
                sum += 2 * deviant;
            }
        }
    }
    *chi_square_ = REAL(sum);
}

void LMFitCPP::calc_model()
{
	std::vector<REAL> & curve = curve_;
	std::vector<REAL> & derivatives = derivatives_;

	calc_curve_values(curve, derivatives);
}
    
void LMFitCPP::calc_coefficients()
{
    std::vector<REAL> & curve = curve_;
    std::vector<REAL> & derivatives = derivatives_;

    calc_chi_square(curve);

    if ((*chi_square_) < prev_chi_square_ || prev_chi_square_ == 0)
    {
        calculate_hessian(derivatives, curve);
        calc_gradient(derivatives, curve);
    }
}

void LMFitCPP::solve_equation_system_gj()
{
    delta_ = gradient_;

    std::vector<REAL> & alpha = modified_hessian_;
    std::vector<REAL> & beta = delta_;

    int icol, irow;
    REAL big, dum, pivinv;

    std::vector<int> indxc(info_.n_parameters_to_fit_, 0);
    std::vector<int> indxr(info_.n_parameters_to_fit_, 0);
    std::vector<int> ipiv(info_.n_parameters_to_fit_, 0);

    for (int kp = 0; kp < info_.n_parameters_to_fit_; kp++)
    {
        big = 0.0;
        for (int jp = 0; jp < info_.n_parameters_to_fit_; jp++)
        {
            if (ipiv[jp] != 1)
            {
                for (int ip = 0; ip < info_.n_parameters_to_fit_; ip++)
                {
                    if (ipiv[ip] == 0)
                    {
                        if (fabs(alpha[jp*info_.n_parameters_to_fit_ + ip]) >= big)
                        {
                            big = fabs(alpha[jp*info_.n_parameters_to_fit_ + ip]);
                            irow = jp;
                            icol = ip;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);


        if (irow != icol)
        {
            for (int ip = 0; ip < info_.n_parameters_to_fit_; ip++)
            {
                std::swap(alpha[irow*info_.n_parameters_to_fit_ + ip], alpha[icol*info_.n_parameters_to_fit_ + ip]);
            }
            std::swap(beta[irow], beta[icol]);
        }
        indxr[kp] = irow;
        indxc[kp] = icol;
        if (alpha[icol*info_.n_parameters_to_fit_ + icol] == 0)
        {
            *state_ = FitState::SINGULAR_HESSIAN;
            break;
        }
        pivinv = 1.f / alpha[icol*info_.n_parameters_to_fit_ + icol];
        alpha[icol*info_.n_parameters_to_fit_ + icol] = 1.0;
        for (int ip = 0; ip < info_.n_parameters_to_fit_; ip++)
        {
            alpha[icol*info_.n_parameters_to_fit_ + ip] *= pivinv;
        }
        beta[icol] *= pivinv;

        for (int jp = 0; jp < info_.n_parameters_to_fit_; jp++)
        {
            if (jp != icol)
            {
                dum = alpha[jp*info_.n_parameters_to_fit_ + icol];
                alpha[jp*info_.n_parameters_to_fit_ + icol] = 0;
                for (int ip = 0; ip < info_.n_parameters_to_fit_; ip++)
                {
                    alpha[jp*info_.n_parameters_to_fit_ + ip] -= alpha[icol*info_.n_parameters_to_fit_ + ip] * dum;
                }
                beta[jp] -= beta[icol] * dum;
            }
        }
    }
}

void LMFitCPP::solve_equation_system_lup()
{
    decompose_hessian_LUP(modified_hessian_);

    solve_LUP(decomposed_hessian_, pivot_array_, gradient_, info_.n_parameters_to_fit_, delta_);
}

void LMFitCPP::update_parameters()
{
    for (int parameter_index = 0, delta_index = 0; parameter_index < info_.n_parameters_; parameter_index++)
    {
        if (parameters_to_fit_[parameter_index])
        {
            prev_parameters_[parameter_index] = parameters_[parameter_index];
            parameters_[parameter_index] = parameters_[parameter_index] + delta_[delta_index++];
        }
    }
}

bool LMFitCPP::check_for_convergence()
{
    bool const fit_found
        = std::abs(*chi_square_ - prev_chi_square_)  < std::max(tolerance_, tolerance_ * std::abs(*chi_square_));

    return fit_found;
}

void LMFitCPP::evaluate_iteration(int const iteration)
{
    bool const max_iterations_reached = iteration == info_.max_n_iterations_ - 1;
    if (converged_ || max_iterations_reached)
    {
        (*n_iterations_) = iteration + 1;
        if (!converged_)
        {
            *state_ = FitState::MAX_ITERATION;
        }
    }
}

void LMFitCPP::prepare_next_iteration()
{
    if ((*chi_square_) < prev_chi_square_)
    {
        lambda_ *= 0.1f;
        prev_chi_square_ = (*chi_square_);
    }
    else
    {
        lambda_ *= 10.;
        (*chi_square_) = prev_chi_square_;
        for (int parameter_index = 0, delta_index = 0; parameter_index < info_.n_parameters_; parameter_index++)
        {
            if (parameters_to_fit_[parameter_index])
            {
                parameters_[parameter_index] = prev_parameters_[parameter_index];
            }
        }
    }
}

void LMFitCPP::modify_step_width()
{
    modified_hessian_ = hessian_;
    size_t const n_parameters = (size_t)(sqrt((REAL)(hessian_.size())));
    for (size_t parameter_index = 0; parameter_index < n_parameters; parameter_index++)
    {
        size_t const diagonal_index = parameter_index * n_parameters + parameter_index;

        // adaptive scaling
        scaling_vector_[parameter_index]
            = std::max(scaling_vector_[parameter_index], modified_hessian_[diagonal_index]);

        // continuous scaling
        //scaling_vector_[parameter_index] = modified_hessian_[diagonal_index];

        // initial scaling
        //if (scaling_vector_[parameter_index] == 0.)
        //    scaling_vector_[parameter_index] = modified_hessian_[diagonal_index];

        modified_hessian_[diagonal_index] += scaling_vector_[parameter_index] * lambda_;
    }
}

void LMFitCPP::run()
{
    for (int i = 0; i < info_.n_parameters_; i++)
        parameters_[i] = initial_parameters_[i];

    *state_ = FitState::CONVERGED;
	calc_model();
    calc_coefficients();

    if (info_.n_parameters_to_fit_ == 0)
        return;

    prev_chi_square_ = (*chi_square_);
        
    for (int iteration = 0; (*state_) == 0; iteration++)
    {
        modify_step_width();
        
        SOLVE_EQUATION_SYSTEM();

        update_parameters();

		calc_model();
        calc_coefficients();

        converged_ = check_for_convergence();

        evaluate_iteration(iteration);

        prepare_next_iteration();

        if (converged_ || *state_ != FitState::CONVERGED)
        {
            break;
        }
    }
}
