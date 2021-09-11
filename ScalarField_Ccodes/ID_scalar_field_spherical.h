
    void ID_scalarfield_spherical(
                     const REAL xyz_or_rthph[3], 
                     const ID_inputs other_inputs,
                     REAL *sf, REAL *sfM) {

      const REAL r  = xyz_or_rthph[0];

      REAL sf_star,psi4_star,alpha_star;

      scalar_field_interpolate_1D(r,
                                  other_inputs.interp_stencil_size,  
                                  other_inputs.numlines_in_file,
                                  other_inputs.r_arr,
                                  other_inputs.sf_arr,
                                  other_inputs.psi4_arr,
                                  other_inputs.alpha_arr,
                                  &sf_star,&psi4_star,&alpha_star);

      // Update varphi
      *sf  = sf_star;
      // Update Pi
      *sfM = 0;

}
