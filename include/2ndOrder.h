// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/* Contains the coefficients of the second order diagonal norm derivative operator of
@article{mattsson2004summation,
  title={Summation by parts operators for finite difference approximations of
         second derivatives},
  author={Mattsson, Ken and Nordstr{\"o}m, Jan},
  journal={Journal of Computational Physics},
  volume={199},
  number={2},
  pages={503--540},
  year={2004},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2004.03.001}
}
*/

#ifdef USE_PERIODIC
//------------------------------------------------------------------------------
// Periodic Derivative Operator
//------------------------------------------------------------------------------
#define ORDER 2
#define NUM_BOUNDS 0
#define STENCIL_WIDTH 3

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // central coefficients
  {-1.0/2.0, 0.0, 1.0/2.0},
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // central coefficients
  1.0,
};

#else
//------------------------------------------------------------------------------
// Non-Periodic Derivative Operator
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS 1
#define STENCIL_WIDTH 3

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, -1.0, 1.0},
  // central coefficients
  {-1.0/2.0, 0.0, 1.0/2.0},
  // right boundary coefficients
  {-1.0, 1.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  2.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  2.0
};

#endif // USE_PERIODIC
