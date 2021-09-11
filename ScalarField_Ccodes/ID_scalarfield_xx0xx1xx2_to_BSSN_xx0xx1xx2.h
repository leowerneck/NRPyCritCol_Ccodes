void ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2(const REAL xx0xx1xx2[3],ID_inputs other_inputs,
                    REAL *sf, REAL *sfM ) {

      REAL sfSphorCart,sfMSphorCart;
      const REAL xx0 = xx0xx1xx2[0];
      const REAL xx1 = xx0xx1xx2[1];
      const REAL xx2 = xx0xx1xx2[2];
      REAL xyz_or_rthph[3];
      xyz_or_rthph[0] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1.0/SINHW) - exp(-1/SINHW));
      xyz_or_rthph[1] = xx1;
      xyz_or_rthph[2] = xx2;
ID_scalarfield_spherical(xyz_or_rthph, other_inputs,
                      &sfSphorCart, &sfMSphorCart);
        // Next compute all rescaled BSSN curvilinear quantities:
      *sf = sfSphorCart;
      *sfM = sfMSphorCart;
}
