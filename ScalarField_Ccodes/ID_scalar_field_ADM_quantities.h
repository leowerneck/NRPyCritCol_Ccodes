
// This function takes as input either (x,y,z) or (r,th,ph) and outputs
//   all ADM quantities in the Cartesian or Spherical basis, respectively.
void ID_scalar_field_ADM_quantities(
                     const REAL xyz_or_rthph[3], 

                     const ID_inputs other_inputs,

                     REAL *gammaDD00,REAL *gammaDD01,REAL *gammaDD02,REAL *gammaDD11,REAL *gammaDD12,REAL *gammaDD22,
                     REAL *KDD00,REAL *KDD01,REAL *KDD02,REAL *KDD11,REAL *KDD12,REAL *KDD22,
                     REAL *alpha,
                     REAL *betaU0,REAL *betaU1,REAL *betaU2,
                     REAL *BU0,REAL *BU1,REAL *BU2) {

      const REAL r  = xyz_or_rthph[0];
      const REAL th = xyz_or_rthph[1];

      REAL sf_star,psi4_star,alpha_star;

      scalar_field_interpolate_1D(r,
                                  other_inputs.interp_stencil_size,  
                                  other_inputs.numlines_in_file,
                                  other_inputs.r_arr,
                                  other_inputs.sf_arr,
                                  other_inputs.psi4_arr,
                                  other_inputs.alpha_arr,
                                  &sf_star,&psi4_star,&alpha_star);

      // Update alpha
      *alpha = alpha_star;
      // \gamma_{rr} = psi^4
      *gammaDD00 = psi4_star;
      // \gamma_{thth} = psi^4 r^2
      *gammaDD11 = psi4_star*r*r;
      // \gamma_{phph} = psi^4 r^2 sin^2(th)
      *gammaDD22 = psi4_star*r*r*sin(th)*sin(th);

      // All other quantities ARE ZERO:
      *gammaDD01 = 0.0; *gammaDD02 = 0.0;
      /**/              *gammaDD12 = 0.0;

      *KDD00 = 0.0; *KDD01 = 0.0; *KDD02 = 0.0;
      /**/          *KDD11 = 0.0; *KDD12 = 0.0;
      /**/                        *KDD22 = 0.0;

      *betaU0 = 0.0; *betaU1 = 0.0; *betaU2 = 0.0;

      *BU0 = 0.0; *BU1 = 0.0; *BU2 = 0.0;
}
