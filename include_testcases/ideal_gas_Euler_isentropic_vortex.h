// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//--------------------------------------------------------------------------------------------------
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
An isentropic vortex, cf. pp. 60f. of
@techreport{shu1997essentially,
  title={Essentially Non-Oscillatory and Weighted Essentially Non-Oscillatory
         Schemes for Hyperbolic Conservation Laws},
  author={Shu, Chi-Wang},
  institution= {NASA},
  month={11},
  year={1997},
  type={Final Report},
  number={NASA/CR-97-206253},
  address={Institute for Computer Applications in Science and Engineering,
           NASA Langley Research Center, Hampton VA United States}
}
*/

// background flow
#define RHO_INF (REAL)(1)
#define P_INF (REAL)(1)
#define VX_INF (REAL)(1)
#define VY_INF (REAL)(1)
#define VZ_INF (REAL)(0)
#define EPSILON (REAL)(5)
#define ALPHA (REAL)(2.5)


inline REAL T_analytical(uint ix, uint iy, uint iz, REAL time) {

  // Coordinates of (ix,iy,iz)
  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  // The solution is u(t,x) = u_0(x - v_inf*t).
  REAL dx = -time*VX_INF;
  REAL dy = -time*VY_INF;
  REAL dz = -time*VZ_INF;

  // A periodic solution is desired.
  // Note: fmod(x, y) seems to return a value in [-y,y] (depending on the sign of x) but we want to have a value in [0,y].
  REAL xx = fmod(x + dx - (REAL)XMIN, (REAL)(XMAX - XMIN)); xx = xx + (xx < 0) * (REAL)(XMAX - XMIN) + (REAL)XMIN;
  REAL yy = fmod(y + dy - (REAL)YMIN, (REAL)(YMAX - YMIN)); yy = yy + (yy < 0) * (REAL)(YMAX - YMIN) + (REAL)YMIN;
  REAL zz = fmod(z + dz - (REAL)ZMIN, (REAL)(ZMAX - ZMIN)); zz = zz + (zz < 0) * (REAL)(ZMAX - ZMIN) + (REAL)ZMIN;

  // place the initial vortex at the centre of the domain
  xx = xx - 0.5 * (XMAX - XMIN);
  yy = yy - 0.5 * (YMAX - YMIN);
  zz = zz - 0.5 * (ZMAX - ZMIN);

  REAL r2 = ALPHA * (xx*xx + yy*yy);

  return P_INF/RHO_INF - (GAMMA-1)/GAMMA * EPSILON*EPSILON / (8*M_PI*M_PI) * exp(-r2);
}


/*
Analytical solution of the density.
*/
inline REAL rho_analytical(uint ix, uint iy, uint iz, REAL time) {

	return RHO_INF * pow(T_analytical(ix, iy, iz, time) * RHO_INF / P_INF, (REAL)(1) / (GAMMA - 1));
}

/*
Analytical solution of the velocity.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

  // Coordinates of (ix,iy,iz)
  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  // The solution is u(t,x) = u_0(x - v_inf*t).
  REAL dx = -time*VX_INF - M_PI;
  REAL dy = -time*VY_INF - M_PI;
  REAL dz = -time*VZ_INF - M_PI;

  // A periodic solution is desired.
  // Note: fmod(x, y) seems to return a value in [-y,y] (depending on the sign of x) but we want to have a value in [0,y].
  REAL xx = fmod(x + dx - (REAL)XMIN, (REAL)(XMAX - XMIN)); xx = xx + (xx < 0) * (REAL)(XMAX - XMIN) + (REAL)XMIN;
  REAL yy = fmod(y + dy - (REAL)YMIN, (REAL)(YMAX - YMIN)); yy = yy + (yy < 0) * (REAL)(YMAX - YMIN) + (REAL)YMIN;
  REAL zz = fmod(z + dz - (REAL)ZMIN, (REAL)(ZMAX - ZMIN)); zz = zz + (zz < 0) * (REAL)(ZMAX - ZMIN) + (REAL)ZMIN;

  // place the initial vortex at the centre of the domain
  xx = xx - 0.5 * (XMAX - XMIN);
  yy = yy - 0.5 * (YMAX - YMIN);
  zz = zz - 0.5 * (ZMAX - ZMIN);

  REAL r2 = ALPHA * (xx*xx + yy*yy);

  REAL vx = VX_INF - yy * EPSILON / (2 * M_PI) * exp(-0.5 * r2);
  REAL vy = VY_INF + xx * EPSILON / (2 * M_PI) * exp(-0.5 * r2);
  REAL vz = VZ_INF;

	return (REAL4) {vx, vy, vz, 0};
}

/*
Analytical solution of the pressure.
*/
inline REAL p_analytical(uint ix, uint iy, uint iz, REAL time) {

	return T_analytical(ix, iy, iz, time) * rho_analytical(ix, iy, iz, time);
}


/*
Boundary condition of the density.
*/
inline REAL rho_boundary(uint ix, uint iy, uint iz, REAL time) {

	return rho_analytical(ix, iy, iz, time);
}

/*
Boundary condition of the velocity.
*/
inline REAL4 u_boundary(uint ix, uint iy, uint iz, REAL time) {

	return u_analytical(ix, iy, iz, time);
}

/*
Boundary condition of the pressure.
*/
inline REAL p_boundary(uint ix, uint iy, uint iz, REAL time) {

	return p_analytical(ix, iy, iz, time);
}


/*
Initial condition of the density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return rho_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the pressure.
*/
inline REAL p_init(uint ix, uint iy, uint iz) {

	return p_analytical(ix, iy, iz, (REAL)(0));
}
