#ifndef GPUFIT_CONSTANTS_H_INCLUDED
#define GPUFIT_CONSTANTS_H_INCLUDED

// fitting model ID

enum ModelID {
    GAUSS_1D = 0,
    GAUSS_2D = 1,
    GAUSS_2D_ELLIPTIC = 2,
    GAUSS_2D_ROTATED = 3,
    CAUCHY_2D_ELLIPTIC = 4,
    LINEAR_1D = 5,
    FLETCHER_POWELL_HELIX = 6,
    BROWN_DENNIS = 7,
    EXP_1D = 8,
    EXP_BASE_1D = 9,
    ERFC_XGAUSS = 10
};

// estimator ID
enum EstimatorID { LSE = 0, MLE = 1 };

// fit state
enum FitState { CONVERGED = 0, MAX_ITERATION = 1, SINGULAR_HESSIAN = 2, NEG_CURVATURE_MLE = 3, GPU_NOT_READY = 4 };

// return state
enum ReturnState { OK = 0, ERROR = -1 };

enum DataLocation { HOST = 0, DEVICE = 1 };

#endif
