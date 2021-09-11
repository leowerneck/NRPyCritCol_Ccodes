void ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(const REAL xx0xx1xx2[3],ID_inputs other_inputs,
                    REAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,
                    REAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,
                    REAL *trK, 
                    REAL *vetU0,REAL *vetU1,REAL *vetU2,
                    REAL *betU0,REAL *betU1,REAL *betU2,
                    REAL *alpha,  REAL *cf) {
      REAL gammaSphorCartDD00,gammaSphorCartDD01,gammaSphorCartDD02,
           gammaSphorCartDD11,gammaSphorCartDD12,gammaSphorCartDD22;
      REAL KSphorCartDD00,KSphorCartDD01,KSphorCartDD02,
           KSphorCartDD11,KSphorCartDD12,KSphorCartDD22;
      REAL alphaSphorCart,betaSphorCartU0,betaSphorCartU1,betaSphorCartU2;
      REAL BSphorCartU0,BSphorCartU1,BSphorCartU2;
      const REAL xx0 = xx0xx1xx2[0];
      const REAL xx1 = xx0xx1xx2[1];
      const REAL xx2 = xx0xx1xx2[2];
      REAL xyz_or_rthph[3];
      xyz_or_rthph[0] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1.0/SINHW) - exp(-1/SINHW));
      xyz_or_rthph[1] = xx1;
      xyz_or_rthph[2] = xx2;
      ID_scalar_field_ADM_quantities(xyz_or_rthph, other_inputs,
                      &gammaSphorCartDD00,&gammaSphorCartDD01,&gammaSphorCartDD02,
                      &gammaSphorCartDD11,&gammaSphorCartDD12,&gammaSphorCartDD22,
                      &KSphorCartDD00,&KSphorCartDD01,&KSphorCartDD02,
                      &KSphorCartDD11,&KSphorCartDD12,&KSphorCartDD22,
                      &alphaSphorCart,&betaSphorCartU0,&betaSphorCartU1,&betaSphorCartU2,
                      &BSphorCartU0,&BSphorCartU1,&BSphorCartU2);
        // Next compute all rescaled BSSN curvilinear quantities:
      const double tmp0 = 1.0/SINHW;
      const double tmp1 = tmp0*xx0;
      const double tmp2 = exp(tmp1);
      const double tmp3 = exp(-tmp1);
      const double tmp4 = tmp0*tmp2 + tmp0*tmp3;
      const double tmp5 = pow(tmp4, 2);
      const double tmp6 = 1.0/tmp5;
      const double tmp7 = pow(AMPL, 2);
      const double tmp8 = exp(tmp0) - exp(-tmp0);
      const double tmp9 = pow(tmp8, 2);
      const double tmp10 = 1.0/tmp9;
      const double tmp11 = tmp10*tmp7;
      const double tmp12 = tmp11*tmp5;
      const double tmp13 = gammaSphorCartDD11*gammaSphorCartDD22;
      const double tmp14 = gammaSphorCartDD00*tmp12;
      const double tmp15 = gammaSphorCartDD01*gammaSphorCartDD02*tmp12;
      const double tmp16 = pow(gammaSphorCartDD12, 2);
      const double tmp17 = pow(gammaSphorCartDD01, 2);
      const double tmp18 = tmp12*tmp17;
      const double tmp19 = pow(gammaSphorCartDD02, 2);
      const double tmp20 = tmp12*tmp19;
      const double tmp21 = -gammaSphorCartDD11*tmp20 + 2*gammaSphorCartDD12*tmp15 - gammaSphorCartDD22*tmp18 + tmp13*tmp14 - tmp14*tmp16;
      const double tmp22 = 1.0/tmp21;
      const double tmp23 = tmp2 - tmp3;
      const double tmp24 = tmp22*pow(tmp23, 4)*tmp5/pow(tmp8, 6);
      const double tmp25 = sin(xx1);
      const double tmp26 = cbrt(tmp24)*pow(fabs(tmp25), 2.0/3.0);
      const double tmp27 = tmp9/tmp7;
      const double tmp28 = 1.0/tmp23;
      const double tmp29 = AMPL*tmp26*tmp28*tmp8;
      const double tmp30 = 1.0/tmp25;
      const double tmp31 = pow(tmp23, 2);
      const double tmp32 = tmp11*tmp31;
      const double tmp33 = tmp26*tmp7;
      const double tmp34 = 1.0/tmp31;
      const double tmp35 = tmp27*tmp34;
      const double tmp36 = tmp26*tmp9;
      const double tmp37 = tmp34*tmp36;
      const double tmp38 = tmp30*tmp37;
      const double tmp39 = pow(tmp25, 2);
      const double tmp40 = 1.0/tmp39;
      const double tmp41 = KSphorCartDD00*tmp12;
      const double tmp42 = 2*tmp22;
      const double tmp43 = AMPL/tmp8;
      const double tmp44 = tmp4*tmp43;
      const double tmp45 = gammaSphorCartDD12*tmp44;
      const double tmp46 = gammaSphorCartDD01*tmp44;
      const double tmp47 = KSphorCartDD01*tmp44;
      const double tmp48 = gammaSphorCartDD02*tmp44;
      const double tmp49 = KSphorCartDD02*tmp44;
      const double tmp50 = KSphorCartDD11*tmp22*(gammaSphorCartDD22*tmp14 - tmp20) + KSphorCartDD12*tmp42*(-gammaSphorCartDD12*tmp14 + tmp15) + KSphorCartDD22*tmp22*(gammaSphorCartDD11*tmp14 - tmp18) + tmp22*tmp41*(tmp13 - tmp16) + tmp42*tmp47*(gammaSphorCartDD02*tmp45 - gammaSphorCartDD22*tmp46) + tmp42*tmp49*(gammaSphorCartDD01*tmp45 - gammaSphorCartDD11*tmp48);
      const double tmp51 = (1.0/3.0)*tmp50;
      const double tmp52 = tmp28*tmp36/tmp4;
      const double tmp53 = tmp23*tmp43;
      const double tmp54 = tmp25*tmp53;
      const double tmp55 = pow(AMPL, 8);
      *hDD00 = tmp27*tmp6*(pow(AMPL, 4)*gammaSphorCartDD00*tmp10*tmp26*tmp5 - tmp12);
      *hDD01 = gammaSphorCartDD01*tmp29;
      *hDD02 = gammaSphorCartDD02*tmp29*tmp30;
      *hDD11 = tmp35*(gammaSphorCartDD11*tmp33 - tmp32);
      *hDD12 = gammaSphorCartDD12*tmp38;
      *hDD22 = tmp35*tmp40*(gammaSphorCartDD22*tmp33 - tmp32*tmp39);
      *aDD00 = tmp36*tmp6*(-tmp14*tmp51 + tmp41);
      *aDD01 = tmp52*(-tmp46*tmp51 + tmp47);
      *aDD02 = tmp30*tmp52*(-tmp48*tmp51 + tmp49);
      *aDD11 = tmp37*(KSphorCartDD11 - gammaSphorCartDD11*tmp51);
      *aDD12 = tmp38*(KSphorCartDD12 - gammaSphorCartDD12*tmp51);
      *aDD22 = tmp37*tmp40*(KSphorCartDD22 - gammaSphorCartDD22*tmp51);
      *trK = tmp50;
      *vetU0 = betaSphorCartU0;
      *vetU1 = betaSphorCartU1*tmp53;
      *vetU2 = betaSphorCartU2*tmp54;
      *betU0 = BSphorCartU0;
      *betU1 = BSphorCartU1*tmp53;
      *betU2 = BSphorCartU2*tmp54;
      *alpha = alphaSphorCart;
      *cf = (1.0/12.0)*log(tmp21/(gammaSphorCartDD00*gammaSphorCartDD11*gammaSphorCartDD22*tmp10*tmp24*tmp39*tmp5*tmp55 - gammaSphorCartDD00*tmp10*tmp16*tmp24*tmp39*tmp5*tmp55 + 2*gammaSphorCartDD01*gammaSphorCartDD02*gammaSphorCartDD12*tmp10*tmp24*tmp39*tmp5*tmp55 - gammaSphorCartDD11*tmp10*tmp19*tmp24*tmp39*tmp5*tmp55 - gammaSphorCartDD22*tmp10*tmp17*tmp24*tmp39*tmp5*tmp55));
}
