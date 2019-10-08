#ifndef GPUFIT_EXP1D_CUH_INCLUDED
#define GPUFIT_EXP1D_CUH_INCLUDED

__device__ void calculate_exp1d(REAL const* parameters,
				int const n_fits,
				int const n_points,
				REAL* value,
				REAL* derivative,
				int const point_index,
				int const fit_index,
				int const chunk_index,
				char* user_info,
				std::size_t const user_info_size){

  REAL* user_info_float = (REAL*) user_info;
  REAL x = 0;
  if(!user_info_float) x = point_index;
  else if(user_info_size / sizeof(REAL) == n_points)
    x = user_info_float[point_index];
  else if(user_info_size / sizeof(REAL) > n_points){
    int const chunk_begin = chunk_index * n_fits * n_points;
    int const fit_begin = fit_index * n_points;
    x = user_info_float[chunk_begin + fit_begin + point_index];
  }

  float ev = exp(-x*parameters[2]);
  value[point_index] = parameters[0] + parameters[1] * ev;

  REAL* current_derivatives = derivative + point_index;
  current_derivatives[0 * n_points] = 1;
  current_derivatives[1 * n_points] = ev;
  current_derivatives[2 * n_points] = -x*ev;
}

#endif
