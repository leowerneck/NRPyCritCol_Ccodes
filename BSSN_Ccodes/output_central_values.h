void output_central_values( const REAL t, const int Nxx_plus_2NGHOSTS[3], REAL *in_gfs, REAL *aux_gfs) {
      
    /* Set indices */
    const int i0 = NGHOSTS;
    const int i1 = NGHOSTS;
    const int i2 = NGHOSTS;
    /* Set needed values of the scalar field */
    const REAL sf_i0p1     = in_gfs[IDX4(SFGF, i0,i1,i2)];
    const REAL sf_i0p2     = in_gfs[IDX4(SFGF, i0+1,i1,i2)];
    const REAL sf_i0p3     = in_gfs[IDX4(SFGF, i0+2,i1,i2)];
    /* Set needed values of alpha */
    const REAL alpha_i0p1  = in_gfs[IDX4(ALPHAGF, i0,i1,i2)];
    const REAL alpha_i0p2  = in_gfs[IDX4(ALPHAGF, i0+1,i1,i2)];
    const REAL alpha_i0p3  = in_gfs[IDX4(ALPHAGF, i0+2,i1,i2)];
    /* Set needed values of T^{\mu\nu} */
    const REAL T4UU00_i0p1 = aux_gfs[IDX4(T4UU00GF, i0,i1,i2)];
    const REAL T4UU00_i0p2 = aux_gfs[IDX4(T4UU00GF, i0+1,i1,i2)];
    const REAL T4UU00_i0p3 = aux_gfs[IDX4(T4UU00GF, i0+2,i1,i2)];
    /* Compute needed values of rho, the energy density */
    const REAL rho_i0p1    = pow(alpha_i0p1,2) * T4UU00_i0p1;
    const REAL rho_i0p2    = pow(alpha_i0p2,2) * T4UU00_i0p2;
    const REAL rho_i0p3    = pow(alpha_i0p3,2) * T4UU00_i0p3;
    /* Compute the central values of the scalar field, alpha, and rho */
    const REAL sf_c    = 3.0*sf_i0p1    - 3.0*sf_i0p2    + sf_i0p3;
    const REAL alpha_c = 3.0*alpha_i0p1 - 3.0*alpha_i0p2 + alpha_i0p3;
    const REAL rho_c   = 3.0*rho_i0p1   - 3.0*rho_i0p2   + rho_i0p3;

    /* Set the output file */
    FILE *outfile;
    if( t > 0.0 ) {
      outfile = fopen("out_central_values.dat","a");
    }
    else {
      outfile = fopen("out_central_values.dat","w");
    }

    /* Output the central values of the scalar field, alpha, and rho */
    fprintf(outfile,"%.15e %.15e %.15e %.15e\n",t,sf_c,alpha_c,rho_c);

    /* Close the file */
    fclose(outfile);

}
