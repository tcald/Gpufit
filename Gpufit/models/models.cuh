#ifndef GPUFIT_MODELS_CUH_INCLUDED
#define GPUFIT_MODELS_CUH_INCLUDED

#include "linear_1d.cuh"
#include "exp_1d.cuh"
#include "exp_base_1d.cuh"
#include "erfc_xgauss.cuh"
#include "gauss_1d.cuh"
#include "gauss_2d.cuh"
#include "gauss_2d_elliptic.cuh"
#include "gauss_2d_rotated.cuh"
#include "cauchy_2d_elliptic.cuh"
#include "fletcher_powell_helix.cuh"
#include "brown_dennis.cuh"

__device__ void calculate_model(
    ModelID const model_id,
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    int const user_info_size)
{
    switch (model_id)
    {
    case GAUSS_1D:
        calculate_gauss1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case GAUSS_2D:
        calculate_gauss2d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case GAUSS_2D_ELLIPTIC:
        calculate_gauss2delliptic(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case GAUSS_2D_ROTATED:
        calculate_gauss2drotated(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case CAUCHY_2D_ELLIPTIC:
        calculate_cauchy2delliptic(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case LINEAR_1D:
        calculate_linear1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case EXP_1D:
      calculate_exp1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
      break;
    case EXP_BASE_1D:
      calculate_expbase1d(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
      break;
    case ERFC_XGAUSS:
      calculate_erfc_xgauss(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
      break;
    case FLETCHER_POWELL_HELIX:
        calculate_fletcher_powell_helix(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    case BROWN_DENNIS:
        calculate_brown_dennis(parameters, n_fits, n_points, value, derivative, point_index, fit_index, chunk_index, user_info, user_info_size);
        break;
    default:
        break;
    }
}

void configure_model(ModelID const model_id, int & n_parameters, int & n_dimensions)
{
    switch (model_id)
    {
    case GAUSS_1D:              n_parameters = 4; n_dimensions = 1; break;
    case GAUSS_2D:              n_parameters = 5; n_dimensions = 2; break;
    case GAUSS_2D_ELLIPTIC:     n_parameters = 6; n_dimensions = 2; break;
    case GAUSS_2D_ROTATED:      n_parameters = 7; n_dimensions = 2; break;
    case CAUCHY_2D_ELLIPTIC:    n_parameters = 6; n_dimensions = 2; break;
    case LINEAR_1D:             n_parameters = 2; n_dimensions = 1; break;
    case EXP_1D:                n_parameters = 2; n_dimensions = 1; break;
    case EXP_BASE_1D:           n_parameters = 3; n_dimensions = 1; break;
    case ERFC_XGAUSS:           n_parameters = 5; n_dimensions = 1; break;
    case FLETCHER_POWELL_HELIX:       n_parameters = 3; n_dimensions = 1; break;
    case BROWN_DENNIS:          n_parameters = 4; n_dimensions = 1; break;
    default:                                                        break;
    }
}

#endif // GPUFIT_MODELS_CUH_INCLUDED
