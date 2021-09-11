int check_regrid_criterion( const int Nxx_plus_2NGHOSTS[3], const REAL dxx[3], REAL *xx[3], REAL *in_gfs ) {

  const REAL inv_dx0      = 1.0/dxx[0];
  const REAL inv_dx0_sqrd = inv_dx0 * inv_dx0;
  const REAL inv_sinhW    = 1.0/SINHW;
  const REAL dr_min       = AMPL * sinh( dxx[0] * inv_sinhW ) / sinh( 1.0 * inv_sinhW );

  for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {

    const REAL phi_i0m2 = in_gfs[IDX4(SFGF, i0-2,NGHOSTS,NGHOSTS)];
    const REAL phi_i0m1 = in_gfs[IDX4(SFGF, i0-1,NGHOSTS,NGHOSTS)];
    const REAL phi      = in_gfs[IDX4(SFGF, i0  ,NGHOSTS,NGHOSTS)];
    const REAL phi_i0p1 = in_gfs[IDX4(SFGF, i0+1,NGHOSTS,NGHOSTS)];
    const REAL phi_i0p2 = in_gfs[IDX4(SFGF, i0+2,NGHOSTS,NGHOSTS)];
    /* First, compute phi_{x} and phi_{xx} */
    const REAL phi_x  = inv_dx0 * ( (1.0/12.0)*( phi_i0m2 - phi_i0p2 ) + (2.0/3.0)*( phi_i0p1 - phi_i0m1 ) );
    const REAL phi_xx = inv_dx0_sqrd * ( -(1.0/12.0)*( phi_i0m2 + phi_i0p2 ) + (4.0/3.0)*( phi_i0p1 + phi_i0m1 ) - (5.0/2.0)*phi );
    /* Then, compute phi_{rr}. We know that
     *
     * r = A sinh( x/w ) / sinh( 1/w )
     *
     * which implies
     *
     * => dr = [ (A/w) cosh( x/w ) / sinh( 1/w ) ] dx
     *
     * => pd_{r} = [ (w/A) sinh( 1/w ) / cosh( x/w ) ] pd_{x}
     *
     * => pd_{rr} = pd_{r} ( pd_{r} )
     *
     * => pd_{rr} = ( (w/A) sinh(1/w) )^{2} ( 1/cosh(x/w) ) pd_{x} [ ( 1/cosh(x/w) ) pd_{x} )
     *
     * => pd_{rr} = ( (w/A) sinh(1/w) )^{2} ( 1/cosh(x/w) )^{2} [ pd_{xx} - tanh(x/w)/w pd_{x} )
     */
    const REAL tmp0   = cosh( xx[0][i0] * inv_sinhW );
    const REAL tmp1   = tanh( xx[0][i0] * inv_sinhW ) * inv_sinhW;
    const REAL tmp2   = 1.0 / pow(tmp0,2.0);
    const REAL tmp3   = pow( SINHW/AMPL * sinh( inv_sinhW ),2.0 );
    const REAL phi_rr = tmp3 * tmp2 * ( phi_xx - tmp1 * phi_x );
    const REAL l      = 1.0 / sqrt( abs(phi_rr) );

    /* If the criterion is met, return true */
    if( l < 15*dr_min ) return 1;

  }

  /* If we reach this point, then the criterion has not been met */
  return 0;

}

void regrid_AMPL( const int interp_stencil_size, const REAL AMPL_NEW, const int Nxx_plus_2NGHOSTS[3], REAL *xx[3], REAL *helper_gfs, REAL *in_and_out_gfs ) {

  /* Useful auxiliary variables */
  const REAL inv_SINHW          = 1.0/SINHW;
  const REAL sinh_inv_SINHW     = sinh( inv_SINHW );
  const REAL inv_sinh_inv_SINHW = 1.0 / sinh_inv_SINHW;

  for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {

    /* Find the new value of r after the regrid */
    const REAL r_star   = AMPL_NEW * sinh( xx[0][i0] * inv_SINHW ) * inv_sinh_inv_SINHW;

    /* Map the new value of r onto the old x grid */
    const REAL x_star   = SINHW * asinh( r_star * sinh_inv_SINHW / AMPL );

    /* Find the index in the x[0] array such that | x[0][idx] - x_star | is minimal */
    const int idx = bisection_idx_finder(x_star, Nxx_plus_2NGHOSTS[0], xx[0]);

    /* Set up minimum and maximum interpolation indices */
    const int idxmin = MAX(0,idx-interp_stencil_size/2-1);
    const int idxmax = idxmin + interp_stencil_size;

    /* printf("x_star: %e | x[0][%d]: %e\n",x_star,idxmin,x[0][idxmin]); */

    /* Compute l_i(x) for the Lagrange polynomial interpolation */
    REAL l_i_of_x[interp_stencil_size];
    for(int i=idxmin;i<idxmax;i++) {
      REAL numer = 1.0;
      REAL denom = 1.0;
      for(int j=idxmin;j<i;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      for(int j=i+1;j<idxmax;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      l_i_of_x[i-idxmin] = numer/denom;      
    }

    /* Perform the interpolation */
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
      for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
	for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
	  const int interp_idx   = IDX4(which_gf,i0,i1,i2);
	  helper_gfs[interp_idx] = 0.0;
	  for( int i=idxmin; i<idxmax; i++ ) {
	    helper_gfs[interp_idx] += l_i_of_x[i-idxmin] * in_and_out_gfs[IDX4(which_gf,i,i1,i2)];
	  }
	}
      }
    }

  } // END OF for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++)

  /* Now update the gridfunctions */
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
      for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
	for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {
	  const int index       = IDX4(which_gf,i0,i1,i2);
	  in_and_out_gfs[index] = helper_gfs[index];
	}
      }
    }
  }

  /* Update the amplitude */
  AMPL = AMPL_NEW;

}

void regrid_SINHW( const int interp_stencil_size, const REAL SINHW_NEW, const int Nxx_plus_2NGHOSTS[3], REAL *xx[3], REAL *helper_gfs, REAL *in_and_out_gfs ) {

  /* Useful auxiliary variables */
  const REAL inv_SINHW          = 1.0/SINHW;
  const REAL sinh_inv_SINHW     = sinh( inv_SINHW );

  const REAL inv_SINHW_NEW          = 1.0/SINHW_NEW;
  const REAL sinh_inv_SINHW_NEW     = sinh( inv_SINHW_NEW );
  const REAL inv_sinh_inv_SINHW_NEW = 1.0 / sinh_inv_SINHW_NEW;

  for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {

    /* Find the new value of r after the regrid */
    const REAL r_star   = AMPL * sinh( xx[0][i0] * inv_SINHW_NEW ) * inv_sinh_inv_SINHW_NEW;

    /* Map the new value of r onto the old x grid */
    const REAL x_star   = SINHW * asinh( r_star * sinh_inv_SINHW / AMPL );

    /* Find the index in the x[0] array such that | x[0][idx] - x_star | is minimal */
    const int idx = bisection_idx_finder(x_star, Nxx_plus_2NGHOSTS[0], xx[0]);

    /* Set up minimum and maximum interpolation indices */
    const int idxmin = MAX(0,idx-interp_stencil_size/2-1);
    const int idxmax = idxmin + interp_stencil_size;

    /* printf("x_star: %e | x[0][%d]: %e\n",x_star,idxmin,x[0][idxmin]); */

    /* Compute l_i(x) for the Lagrange polynomial interpolation */
    REAL l_i_of_x[interp_stencil_size];
    for(int i=idxmin;i<idxmax;i++) {
      REAL numer = 1.0;
      REAL denom = 1.0;
      for(int j=idxmin;j<i;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      for(int j=i+1;j<idxmax;j++) {
	numer *= x_star   - xx[0][j];
	denom *= xx[0][i] - xx[0][j];
      }
      l_i_of_x[i-idxmin] = numer/denom;      
    }

    /* Perform the interpolation */
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
      for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
	for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
	  const int interp_idx   = IDX4(which_gf,i0,i1,i2);
	  helper_gfs[interp_idx] = 0.0;
	  for( int i=idxmin; i<idxmax; i++ ) {
	    helper_gfs[interp_idx] += l_i_of_x[i-idxmin] * in_and_out_gfs[IDX4(which_gf,i,i1,i2)];
	  }
	}
      }
    }

  } // END OF for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++)

  /* Now update the gridfunctions */
  for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
    for(int i2=NGHOSTS;i2<Nxx_plus_2NGHOSTS[2]-NGHOSTS;i2++) {
      for(int i1=NGHOSTS;i1<Nxx_plus_2NGHOSTS[1]-NGHOSTS;i1++) {
	for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS[0]-NGHOSTS;i0++) {
	  const int index       = IDX4(which_gf,i0,i1,i2);
	  in_and_out_gfs[index] = helper_gfs[index];
	}
      }
    }
  }

  /* Update SINHW */
  SINHW = SINHW_NEW;

}
