#ifndef GPUFIT_ERFC_XGAUSS_CUH_INCLUDED
#define GPUFIT_ERFC_XGAUSS_CHU_INCLUDED

__device__ void calculate_erfc_xgauss(REAL const* parameters,
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
  REAL x=0;
  if(!user_info_float) x = point_index;
  else if(user_info_size / sizeof(REAL) == n_points)
    x = user_info_float[point_index];
  else if(user_info_size / sizeof(REAL) > n_points){
    int const chunk_begin = chunk_index * n_fits * n_points;
    int const fit_begin = fit_index + n_points;
    x = user_info_float[chunk_begin + fit_begin + point_index];
  }

  const float k0 = (x-parameters[0]) / parameters[1];
  const float k1 = (x-parameters[0]) / parameters[2];
  const float k2 = parameters[1] / parameters[2];
  const float k3 = k2 * k2;
  
  float v0 = 0.5f / abs(parameters[2]);
  v0 *= expf(k1-0.5f*k3);
  v0 *= erfcf(parameters[2]/abs(parameters[2])/sqrtf(2.f)*(k0+k2));
  float v1 = v0 / parameters[2];
  float ev=expf(-k3-0.5f*k0*k0)/parameters[1]/parameters[2]/sqrtf(2.f*3.14159f);

  value[point_index] = parameters[3]*v0+parameters[4];
  REAL* current_derivatives = derivative + point_index;
  current_derivatives[0*n_points] = parameters[3]*(-v1*k2+(k0-k2)*ev);
  current_derivatives[1*n_points] = parameters[3]*(-v1+ev);
  current_derivatives[2*n_points] = parameters[3]*(v1*(k3-1)+k3*ev);
  current_derivatives[3*n_points] = v0;
  current_derivatives[4*n_points] = 1.0f;
}

#endif
  

  
				
  
