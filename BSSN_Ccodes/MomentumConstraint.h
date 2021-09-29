const REAL invdx0 = 1.0/dxx[0];
for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++) {
    for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++) {
        const REAL xx1 = xx[1][i1];
#pragma omp parallel for
        for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++) {
            const REAL xx0 = xx[0][i0];
            {
               /* 
                * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                */
               const double hDD00_i0m2_i1_i2 = in_gfs[IDX4(HDD00GF, i0-2,i1,i2)];
               const double hDD00_i0m1_i1_i2 = in_gfs[IDX4(HDD00GF, i0-1,i1,i2)];
               const double hDD00 = in_gfs[IDX4(HDD00GF, i0,i1,i2)];
               const double hDD00_i0p1_i1_i2 = in_gfs[IDX4(HDD00GF, i0+1,i1,i2)];
               const double hDD00_i0p2_i1_i2 = in_gfs[IDX4(HDD00GF, i0+2,i1,i2)];
               const double hDD01_i0m2_i1_i2 = in_gfs[IDX4(HDD01GF, i0-2,i1,i2)];
               const double hDD01_i0m1_i1_i2 = in_gfs[IDX4(HDD01GF, i0-1,i1,i2)];
               const double hDD01 = in_gfs[IDX4(HDD01GF, i0,i1,i2)];
               const double hDD01_i0p1_i1_i2 = in_gfs[IDX4(HDD01GF, i0+1,i1,i2)];
               const double hDD01_i0p2_i1_i2 = in_gfs[IDX4(HDD01GF, i0+2,i1,i2)];
               const double hDD02_i0m2_i1_i2 = in_gfs[IDX4(HDD02GF, i0-2,i1,i2)];
               const double hDD02_i0m1_i1_i2 = in_gfs[IDX4(HDD02GF, i0-1,i1,i2)];
               const double hDD02 = in_gfs[IDX4(HDD02GF, i0,i1,i2)];
               const double hDD02_i0p1_i1_i2 = in_gfs[IDX4(HDD02GF, i0+1,i1,i2)];
               const double hDD02_i0p2_i1_i2 = in_gfs[IDX4(HDD02GF, i0+2,i1,i2)];
               const double hDD11_i0m2_i1_i2 = in_gfs[IDX4(HDD11GF, i0-2,i1,i2)];
               const double hDD11_i0m1_i1_i2 = in_gfs[IDX4(HDD11GF, i0-1,i1,i2)];
               const double hDD11 = in_gfs[IDX4(HDD11GF, i0,i1,i2)];
               const double hDD11_i0p1_i1_i2 = in_gfs[IDX4(HDD11GF, i0+1,i1,i2)];
               const double hDD11_i0p2_i1_i2 = in_gfs[IDX4(HDD11GF, i0+2,i1,i2)];
               const double hDD12_i0m2_i1_i2 = in_gfs[IDX4(HDD12GF, i0-2,i1,i2)];
               const double hDD12_i0m1_i1_i2 = in_gfs[IDX4(HDD12GF, i0-1,i1,i2)];
               const double hDD12 = in_gfs[IDX4(HDD12GF, i0,i1,i2)];
               const double hDD12_i0p1_i1_i2 = in_gfs[IDX4(HDD12GF, i0+1,i1,i2)];
               const double hDD12_i0p2_i1_i2 = in_gfs[IDX4(HDD12GF, i0+2,i1,i2)];
               const double hDD22_i0m2_i1_i2 = in_gfs[IDX4(HDD22GF, i0-2,i1,i2)];
               const double hDD22_i0m1_i1_i2 = in_gfs[IDX4(HDD22GF, i0-1,i1,i2)];
               const double hDD22 = in_gfs[IDX4(HDD22GF, i0,i1,i2)];
               const double hDD22_i0p1_i1_i2 = in_gfs[IDX4(HDD22GF, i0+1,i1,i2)];
               const double hDD22_i0p2_i1_i2 = in_gfs[IDX4(HDD22GF, i0+2,i1,i2)];
               const double aDD00_i0m2_i1_i2 = in_gfs[IDX4(ADD00GF, i0-2,i1,i2)];
               const double aDD00_i0m1_i1_i2 = in_gfs[IDX4(ADD00GF, i0-1,i1,i2)];
               const double aDD00 = in_gfs[IDX4(ADD00GF, i0,i1,i2)];
               const double aDD00_i0p1_i1_i2 = in_gfs[IDX4(ADD00GF, i0+1,i1,i2)];
               const double aDD00_i0p2_i1_i2 = in_gfs[IDX4(ADD00GF, i0+2,i1,i2)];
               const double aDD01_i0m2_i1_i2 = in_gfs[IDX4(ADD01GF, i0-2,i1,i2)];
               const double aDD01_i0m1_i1_i2 = in_gfs[IDX4(ADD01GF, i0-1,i1,i2)];
               const double aDD01 = in_gfs[IDX4(ADD01GF, i0,i1,i2)];
               const double aDD01_i0p1_i1_i2 = in_gfs[IDX4(ADD01GF, i0+1,i1,i2)];
               const double aDD01_i0p2_i1_i2 = in_gfs[IDX4(ADD01GF, i0+2,i1,i2)];
               const double aDD02_i0m2_i1_i2 = in_gfs[IDX4(ADD02GF, i0-2,i1,i2)];
               const double aDD02_i0m1_i1_i2 = in_gfs[IDX4(ADD02GF, i0-1,i1,i2)];
               const double aDD02 = in_gfs[IDX4(ADD02GF, i0,i1,i2)];
               const double aDD02_i0p1_i1_i2 = in_gfs[IDX4(ADD02GF, i0+1,i1,i2)];
               const double aDD02_i0p2_i1_i2 = in_gfs[IDX4(ADD02GF, i0+2,i1,i2)];
               const double aDD11_i0m2_i1_i2 = in_gfs[IDX4(ADD11GF, i0-2,i1,i2)];
               const double aDD11_i0m1_i1_i2 = in_gfs[IDX4(ADD11GF, i0-1,i1,i2)];
               const double aDD11 = in_gfs[IDX4(ADD11GF, i0,i1,i2)];
               const double aDD11_i0p1_i1_i2 = in_gfs[IDX4(ADD11GF, i0+1,i1,i2)];
               const double aDD11_i0p2_i1_i2 = in_gfs[IDX4(ADD11GF, i0+2,i1,i2)];
               const double aDD12_i0m2_i1_i2 = in_gfs[IDX4(ADD12GF, i0-2,i1,i2)];
               const double aDD12_i0m1_i1_i2 = in_gfs[IDX4(ADD12GF, i0-1,i1,i2)];
               const double aDD12 = in_gfs[IDX4(ADD12GF, i0,i1,i2)];
               const double aDD12_i0p1_i1_i2 = in_gfs[IDX4(ADD12GF, i0+1,i1,i2)];
               const double aDD12_i0p2_i1_i2 = in_gfs[IDX4(ADD12GF, i0+2,i1,i2)];
               const double aDD22_i0m2_i1_i2 = in_gfs[IDX4(ADD22GF, i0-2,i1,i2)];
               const double aDD22_i0m1_i1_i2 = in_gfs[IDX4(ADD22GF, i0-1,i1,i2)];
               const double aDD22 = in_gfs[IDX4(ADD22GF, i0,i1,i2)];
               const double aDD22_i0p1_i1_i2 = in_gfs[IDX4(ADD22GF, i0+1,i1,i2)];
               const double aDD22_i0p2_i1_i2 = in_gfs[IDX4(ADD22GF, i0+2,i1,i2)];
               const double vetU0 = in_gfs[IDX4(VETU0GF, i0,i1,i2)];
               const double vetU1 = in_gfs[IDX4(VETU1GF, i0,i1,i2)];
               const double vetU2 = in_gfs[IDX4(VETU2GF, i0,i1,i2)];
               const double trK_i0m2_i1_i2 = in_gfs[IDX4(TRKGF, i0-2,i1,i2)];
               const double trK_i0m1_i1_i2 = in_gfs[IDX4(TRKGF, i0-1,i1,i2)];
               const double trK_i0p1_i1_i2 = in_gfs[IDX4(TRKGF, i0+1,i1,i2)];
               const double trK_i0p2_i1_i2 = in_gfs[IDX4(TRKGF, i0+2,i1,i2)];
               const double cf_i0m2_i1_i2 = in_gfs[IDX4(CFGF, i0-2,i1,i2)];
               const double cf_i0m1_i1_i2 = in_gfs[IDX4(CFGF, i0-1,i1,i2)];
               const double cf = in_gfs[IDX4(CFGF, i0,i1,i2)];
               const double cf_i0p1_i1_i2 = in_gfs[IDX4(CFGF, i0+1,i1,i2)];
               const double cf_i0p2_i1_i2 = in_gfs[IDX4(CFGF, i0+2,i1,i2)];
               const double alpha = in_gfs[IDX4(ALPHAGF, i0,i1,i2)];
               const double T4UU00 = aux_gfs[IDX4(T4UU00GF, i0,i1,i2)];
               const double T4UU01 = aux_gfs[IDX4(T4UU01GF, i0,i1,i2)];
               const double T4UU02 = aux_gfs[IDX4(T4UU02GF, i0,i1,i2)];
               const double T4UU03 = aux_gfs[IDX4(T4UU03GF, i0,i1,i2)];
               const double aDD_dD000 = invdx0*(-2.0/3.0*aDD00_i0m1_i1_i2 + (1.0/12.0)*aDD00_i0m2_i1_i2 + (2.0/3.0)*aDD00_i0p1_i1_i2 - 1.0/12.0*aDD00_i0p2_i1_i2);
               const double aDD_dD010 = invdx0*(-2.0/3.0*aDD01_i0m1_i1_i2 + (1.0/12.0)*aDD01_i0m2_i1_i2 + (2.0/3.0)*aDD01_i0p1_i1_i2 - 1.0/12.0*aDD01_i0p2_i1_i2);
               const double aDD_dD020 = invdx0*(-2.0/3.0*aDD02_i0m1_i1_i2 + (1.0/12.0)*aDD02_i0m2_i1_i2 + (2.0/3.0)*aDD02_i0p1_i1_i2 - 1.0/12.0*aDD02_i0p2_i1_i2);
               const double aDD_dD110 = invdx0*(-2.0/3.0*aDD11_i0m1_i1_i2 + (1.0/12.0)*aDD11_i0m2_i1_i2 + (2.0/3.0)*aDD11_i0p1_i1_i2 - 1.0/12.0*aDD11_i0p2_i1_i2);
               const double aDD_dD120 = invdx0*(-2.0/3.0*aDD12_i0m1_i1_i2 + (1.0/12.0)*aDD12_i0m2_i1_i2 + (2.0/3.0)*aDD12_i0p1_i1_i2 - 1.0/12.0*aDD12_i0p2_i1_i2);
               const double aDD_dD220 = invdx0*(-2.0/3.0*aDD22_i0m1_i1_i2 + (1.0/12.0)*aDD22_i0m2_i1_i2 + (2.0/3.0)*aDD22_i0p1_i1_i2 - 1.0/12.0*aDD22_i0p2_i1_i2);
               const double cf_dD0 = invdx0*(-2.0/3.0*cf_i0m1_i1_i2 + (1.0/12.0)*cf_i0m2_i1_i2 + (2.0/3.0)*cf_i0p1_i1_i2 - 1.0/12.0*cf_i0p2_i1_i2);
               const double hDD_dD000 = invdx0*(-2.0/3.0*hDD00_i0m1_i1_i2 + (1.0/12.0)*hDD00_i0m2_i1_i2 + (2.0/3.0)*hDD00_i0p1_i1_i2 - 1.0/12.0*hDD00_i0p2_i1_i2);
               const double hDD_dD010 = invdx0*(-2.0/3.0*hDD01_i0m1_i1_i2 + (1.0/12.0)*hDD01_i0m2_i1_i2 + (2.0/3.0)*hDD01_i0p1_i1_i2 - 1.0/12.0*hDD01_i0p2_i1_i2);
               const double hDD_dD020 = invdx0*(-2.0/3.0*hDD02_i0m1_i1_i2 + (1.0/12.0)*hDD02_i0m2_i1_i2 + (2.0/3.0)*hDD02_i0p1_i1_i2 - 1.0/12.0*hDD02_i0p2_i1_i2);
               const double hDD_dD110 = invdx0*(-2.0/3.0*hDD11_i0m1_i1_i2 + (1.0/12.0)*hDD11_i0m2_i1_i2 + (2.0/3.0)*hDD11_i0p1_i1_i2 - 1.0/12.0*hDD11_i0p2_i1_i2);
               const double hDD_dD120 = invdx0*(-2.0/3.0*hDD12_i0m1_i1_i2 + (1.0/12.0)*hDD12_i0m2_i1_i2 + (2.0/3.0)*hDD12_i0p1_i1_i2 - 1.0/12.0*hDD12_i0p2_i1_i2);
               const double hDD_dD220 = invdx0*(-2.0/3.0*hDD22_i0m1_i1_i2 + (1.0/12.0)*hDD22_i0m2_i1_i2 + (2.0/3.0)*hDD22_i0p1_i1_i2 - 1.0/12.0*hDD22_i0p2_i1_i2);
               const double trK_dD0 = invdx0*(-2.0/3.0*trK_i0m1_i1_i2 + (1.0/12.0)*trK_i0m2_i1_i2 + (2.0/3.0)*trK_i0p1_i1_i2 - 1.0/12.0*trK_i0p2_i1_i2);
               /* 
                * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                */
               const double tmp0 = 1.0/SINHW;
               const double tmp1 = tmp0*xx0;
               const double tmp2 = exp(tmp1);
               const double tmp3 = tmp0*tmp2;
               const double tmp4 = exp(-tmp1);
               const double tmp5 = tmp0*tmp4;
               const double tmp6 = tmp3 + tmp5;
               const double tmp7 = exp(tmp0) - exp(-tmp0);
               const double tmp8 = pow(AMPL, 4)/pow(tmp7, 4);
               const double tmp9 = sin(xx1);
               const double tmp10 = pow(tmp9, 2);
               const double tmp11 = tmp2 - tmp4;
               const double tmp12 = tmp10*pow(tmp11, 4);
               const double tmp13 = pow(hDD12, 2)*tmp12*tmp8;
               const double tmp14 = pow(tmp11, 2);
               const double tmp15 = pow(AMPL, 2);
               const double tmp16 = tmp15/pow(tmp7, 2);
               const double tmp17 = tmp14*tmp16;
               const double tmp18 = hDD11*tmp17 + tmp17;
               const double tmp19 = tmp10*tmp17;
               const double tmp20 = hDD22*tmp19 + tmp19;
               const double tmp21 = tmp18*tmp20;
               const double tmp22 = -tmp13 + tmp21;
               const double tmp23 = pow(tmp6, 2);
               const double tmp24 = 2*hDD01;
               const double tmp25 = pow(AMPL, 6)*hDD02*hDD12*tmp12*tmp23*tmp24/pow(tmp7, 6);
               const double tmp26 = tmp16*tmp23;
               const double tmp27 = hDD00*tmp26 + tmp26;
               const double tmp28 = tmp13*tmp27;
               const double tmp29 = tmp14*tmp23*tmp8;
               const double tmp30 = pow(hDD02, 2)*tmp10*tmp29;
               const double tmp31 = tmp18*tmp30;
               const double tmp32 = pow(hDD01, 2)*tmp29;
               const double tmp33 = tmp20*tmp32;
               const double tmp34 = tmp21*tmp27;
               const double tmp35 = tmp25 - tmp28 - tmp31 - tmp33 + tmp34;
               const double tmp36 = 1.0/tmp35;
               const double tmp37 = (2.0/3.0)*tmp36*trK_dD0;
               const double tmp38 = exp(8*cf);
               const double tmp39 = hDD02*tmp6;
               const double tmp40 = hDD12*pow(tmp11, 3)*tmp8;
               const double tmp41 = tmp10*tmp39*tmp40;
               const double tmp42 = tmp11*tmp16;
               const double tmp43 = tmp42*tmp6;
               const double tmp44 = hDD01*tmp43;
               const double tmp45 = tmp20*tmp44;
               const double tmp46 = tmp38*tmp41 - tmp38*tmp45;
               const double tmp47 = exp(12*cf);
               const double tmp48 = 1.0/(tmp25*tmp47 - tmp28*tmp47 - tmp31*tmp47 - tmp33*tmp47 + tmp34*tmp47);
               const double tmp49 = 4*cf;
               const double tmp50 = exp(tmp49);
               const double tmp51 = alpha*tmp50;
               const double tmp52 = T4UU03*tmp51;
               const double tmp53 = tmp17*tmp9;
               const double tmp54 = hDD12*tmp53;
               const double tmp55 = T4UU01*tmp51;
               const double tmp56 = T4UU02*tmp51;
               const double tmp57 = AMPL/tmp7;
               const double tmp58 = tmp50*tmp57;
               const double tmp59 = tmp11*tmp58*vetU0;
               const double tmp60 = tmp58*vetU2;
               const double tmp61 = hDD12*tmp11;
               const double tmp62 = tmp50*tmp7/AMPL;
               const double tmp63 = tmp62/tmp11;
               const double tmp64 = T4UU00*alpha;
               const double tmp65 = tmp48*(tmp18*tmp56 + tmp44*tmp55 + tmp52*tmp54 - tmp64*(-hDD01*tmp59 - tmp18*tmp63*vetU1 - tmp60*tmp61));
               const double tmp66 = hDD01*tmp6;
               const double tmp67 = tmp40*tmp66*tmp9;
               const double tmp68 = tmp43*tmp9;
               const double tmp69 = hDD02*tmp68;
               const double tmp70 = tmp18*tmp69;
               const double tmp71 = tmp38*tmp67 - tmp38*tmp70;
               const double tmp72 = hDD02*tmp9;
               const double tmp73 = tmp58*vetU1;
               const double tmp74 = tmp48*(tmp20*tmp52 + tmp54*tmp56 + tmp55*tmp69 - tmp64*(-tmp20*tmp63*vetU2/tmp9 - tmp59*tmp72 - tmp61*tmp73*tmp9));
               const double tmp75 = tmp48*(tmp27*tmp55 + tmp44*tmp56 + tmp52*tmp69 - tmp64*(-tmp27*tmp62*vetU0/tmp6 - tmp39*tmp60 - tmp66*tmp73));
               const double tmp76 = 8*M_PI*tmp50;
               const double tmp77 = pow(tmp35, -2);
               const double tmp78 = tmp67 - tmp70;
               const double tmp79 = tmp77*pow(tmp78, 2);
               const double tmp80 = aDD22*tmp19;
               const double tmp81 = 6*tmp80;
               const double tmp82 = aDD11*tmp17;
               const double tmp83 = tmp41 - tmp45;
               const double tmp84 = tmp77*pow(tmp83, 2);
               const double tmp85 = 6*tmp84;
               const double tmp86 = pow(tmp22, 2);
               const double tmp87 = aDD00*tmp26;
               const double tmp88 = 6*tmp87;
               const double tmp89 = tmp77*tmp78;
               const double tmp90 = tmp83*tmp89;
               const double tmp91 = aDD12*tmp53;
               const double tmp92 = aDD02*tmp68;
               const double tmp93 = tmp22*tmp77;
               const double tmp94 = 12*tmp93;
               const double tmp95 = aDD01*tmp43;
               const double tmp96 = hDD01*tmp29*tmp72;
               const double tmp97 = tmp27*tmp54;
               const double tmp98 = tmp96 - tmp97;
               const double tmp99 = cos(xx1);
               const double tmp100 = tmp17*tmp99;
               const double tmp101 = hDD12*tmp100*tmp36;
               const double tmp102 = hDD_dD110*tmp17;
               const double tmp103 = 2*tmp0;
               const double tmp104 = exp(tmp103);
               const double tmp105 = exp(2*tmp1);
               const double tmp106 = tmp15*(tmp105 - 1)*(tmp105 + 1)*exp(-tmp103*(xx0 - 1))/pow(tmp104 - 1, 2);
               const double tmp107 = tmp103*tmp106;
               const double tmp108 = tmp42*(2*tmp3 + 2*tmp5);
               const double tmp109 = hDD11*tmp108;
               const double tmp110 = (1.0/2.0)*tmp36;
               const double tmp111 = tmp110*(-tmp102 - tmp107 - tmp109);
               const double tmp112 = -tmp101*tmp98 - tmp111*tmp83;
               const double tmp113 = -tmp101*tmp78 - tmp111*tmp22;
               const double tmp114 = 2*tmp43;
               const double tmp115 = aDD01*tmp114;
               const double tmp116 = tmp18*tmp27;
               const double tmp117 = tmp116 - tmp32;
               const double tmp118 = -tmp101*tmp117 - tmp111*tmp78;
               const double tmp119 = 2*tmp91;
               const double tmp120 = tmp77*(2*tmp112*tmp82 + tmp113*tmp115 + tmp118*tmp119);
               const double tmp121 = tmp20*tmp27;
               const double tmp122 = tmp121 - tmp30;
               const double tmp123 = tmp122*tmp83;
               const double tmp124 = tmp83*tmp98;
               const double tmp125 = 2*tmp53*tmp99;
               const double tmp126 = hDD22*tmp125;
               const double tmp127 = tmp15*pow(1 - tmp105, 2)*exp(tmp103*(1 - xx0))*sin(2*xx1)/pow(1 - tmp104, 2);
               const double tmp128 = tmp110*(tmp126 + tmp127);
               const double tmp129 = hDD_dD120*tmp53;
               const double tmp130 = tmp43*tmp99;
               const double tmp131 = hDD02*tmp130;
               const double tmp132 = tmp108*tmp9;
               const double tmp133 = hDD12*tmp132;
               const double tmp134 = tmp110*(-tmp129 + tmp131 - tmp133);
               const double tmp135 = -tmp128*tmp98 - tmp134*tmp83;
               const double tmp136 = tmp135*tmp82;
               const double tmp137 = -tmp128*tmp78 - tmp134*tmp22;
               const double tmp138 = -tmp117*tmp128 - tmp134*tmp78;
               const double tmp139 = tmp138*tmp91;
               const double tmp140 = tmp115*tmp137 + 2*tmp136 + 2*tmp139;
               const double tmp141 = tmp140*tmp77;
               const double tmp142 = tmp110*(-tmp126 - tmp127);
               const double tmp143 = hDD_dD220*tmp19;
               const double tmp144 = tmp10*tmp107;
               const double tmp145 = tmp10*tmp108;
               const double tmp146 = hDD22*tmp145;
               const double tmp147 = tmp110*(-tmp143 - tmp144 - tmp146);
               const double tmp148 = -tmp142*tmp98 - tmp147*tmp78;
               const double tmp149 = -tmp142*tmp83 - tmp147*tmp22;
               const double tmp150 = tmp114*tmp9;
               const double tmp151 = aDD02*tmp150;
               const double tmp152 = -tmp122*tmp142 - tmp147*tmp83;
               const double tmp153 = tmp77*(tmp119*tmp152 + 2*tmp148*tmp80 + tmp149*tmp151);
               const double tmp154 = tmp117*tmp153;
               const double tmp155 = tmp135*tmp91;
               const double tmp156 = tmp138*tmp80;
               const double tmp157 = aDD22*tmp125 + tmp137*tmp151 + 2*tmp155 + 2*tmp156;
               const double tmp158 = tmp89*tmp98;
               const double tmp159 = tmp110*(tmp102 + tmp107 + tmp109);
               const double tmp160 = tmp129 + tmp133;
               const double tmp161 = tmp110*(tmp131 + tmp160);
               const double tmp162 = -tmp117*tmp161 - tmp159*tmp98;
               const double tmp163 = tmp162*tmp91;
               const double tmp164 = -tmp159*tmp83 - tmp161*tmp78;
               const double tmp165 = -tmp122*tmp159 - tmp161*tmp98;
               const double tmp166 = tmp165*tmp82;
               const double tmp167 = aDD11*tmp108 + aDD_dD110*tmp17 + tmp115*tmp164 + 2*tmp163 + 2*tmp166;
               const double tmp168 = tmp110*(-tmp131 + tmp160);
               const double tmp169 = tmp110*(tmp143 + tmp144 + tmp146);
               const double tmp170 = -tmp117*tmp169 - tmp168*tmp98;
               const double tmp171 = tmp170*tmp80;
               const double tmp172 = -tmp168*tmp83 - tmp169*tmp78;
               const double tmp173 = -tmp122*tmp168 - tmp169*tmp98;
               const double tmp174 = tmp173*tmp91;
               const double tmp175 = aDD22*tmp145 + aDD_dD220*tmp19 + tmp151*tmp172 + 2*tmp171 + 2*tmp174;
               const double tmp176 = tmp164*tmp87;
               const double tmp177 = tmp115*tmp165 + tmp151*tmp162 + 2*tmp176;
               const double tmp178 = tmp177*tmp83;
               const double tmp179 = tmp172*tmp87;
               const double tmp180 = tmp115*tmp173 + tmp151*tmp170 + 2*tmp179;
               const double tmp181 = tmp180*tmp78;
               const double tmp182 = pow(SINHW, -2);
               const double tmp183 = tmp182*tmp2;
               const double tmp184 = tmp182*tmp4;
               const double tmp185 = tmp16*tmp6*(2*tmp183 - 2*tmp184);
               const double tmp186 = tmp110*(hDD00*tmp185 + hDD_dD000*tmp26 + 2*tmp106/pow(SINHW, 3));
               const double tmp187 = tmp42*(tmp183 - tmp184);
               const double tmp188 = tmp187*tmp9 + tmp26*tmp9;
               const double tmp189 = tmp110*(2*hDD02*tmp188 + hDD_dD020*tmp150);
               const double tmp190 = tmp187 + tmp26;
               const double tmp191 = tmp110*(hDD_dD010*tmp114 + tmp190*tmp24);
               const double tmp192 = -tmp186*tmp22 - tmp189*tmp78 - tmp191*tmp83;
               const double tmp193 = -tmp117*tmp189 - tmp186*tmp78 - tmp191*tmp98;
               const double tmp194 = -tmp122*tmp191 - tmp186*tmp83 - tmp189*tmp98;
               const double tmp195 = tmp77*(aDD00*tmp185 + aDD_dD000*tmp26 + tmp115*tmp194 + tmp151*tmp193 + 2*tmp192*tmp87);
               const double tmp196 = tmp112*tmp95 + tmp113*tmp87 + tmp118*tmp92 + tmp163 + tmp164*tmp95 + tmp166;
               const double tmp197 = aDD12*tmp100 + tmp112*tmp91 + tmp113*tmp92 + tmp118*tmp80 + tmp136 + tmp137*tmp95 + tmp139;
               const double tmp198 = tmp197*tmp77;
               const double tmp199 = tmp122*tmp78;
               const double tmp200 = tmp137*tmp92 + tmp148*tmp91 + tmp149*tmp95 + tmp152*tmp82 + tmp155 + tmp156;
               const double tmp201 = tmp117*tmp77;
               const double tmp202 = tmp200*tmp201;
               const double tmp203 = tmp122*tmp93;
               const double tmp204 = tmp148*tmp92 + tmp149*tmp87 + tmp152*tmp95 + tmp171 + tmp172*tmp92 + tmp174;
               const double tmp205 = tmp135*tmp95 + tmp137*tmp87 + tmp138*tmp92;
               const double tmp206 = tmp170*tmp91 + tmp172*tmp95 + tmp173*tmp82;
               const double tmp207 = tmp205 + tmp206;
               const double tmp208 = tmp93*tmp98;
               const double tmp209 = tmp162*tmp80 + tmp164*tmp92 + tmp165*tmp91;
               const double tmp210 = aDD02*tmp130 + tmp205 + tmp209;
               const double tmp211 = tmp201*tmp204;
               const double tmp212 = aDD12*tmp132 + aDD_dD120*tmp53 + tmp206 + tmp209;
               const double tmp213 = aDD01*tmp190 + aDD_dD010*tmp43 + tmp162*tmp92 + tmp165*tmp95 + tmp176 + tmp192*tmp95 + tmp193*tmp91 + tmp194*tmp82;
               const double tmp214 = 2*tmp93;
               const double tmp215 = aDD02*tmp188 + aDD_dD020*tmp68 + tmp170*tmp92 + tmp173*tmp95 + tmp179 + tmp192*tmp92 + tmp193*tmp80 + tmp194*tmp91;
               const double tmp216 = tmp57*exp(-tmp49);
               const double tmp217 = tmp38*tmp96 - tmp38*tmp97;
               const double tmp218 = tmp77*pow(tmp98, 2);
               const double tmp219 = tmp122*tmp98;
               const double tmp220 = tmp167*tmp77;
               const double tmp221 = tmp124*tmp77;
               const double tmp222 = 6*tmp91;
               const double tmp223 = tmp88*tmp93;
               const double tmp224 = 6*tmp92;
               const double tmp225 = tmp199*tmp77;
               const double tmp226 = tmp123*tmp77;
               const double tmp227 = 6*tmp82;
               const double tmp228 = 6*tmp95;
               const double tmp229 = tmp195*tmp22;
               const double tmp230 = tmp11*tmp216;
               const double tmp231 = tmp201*tmp78;
               const double tmp232 = tmp201*tmp83;
               const double tmp233 = tmp201*tmp22;
               aux_gfs[IDX4(MU0GF, i0, i1, i2)] = tmp216*tmp6*(cf_dD0*(tmp77*tmp86*tmp88 + tmp78*tmp92*tmp94 + tmp79*tmp81 + tmp82*tmp85 + tmp83*tmp94*tmp95 + 12*tmp90*tmp91) + tmp120*tmp123 + tmp124*tmp141 + tmp124*tmp198 + tmp154*tmp78 + tmp157*tmp158 + tmp158*tmp200 + tmp167*tmp84 + tmp175*tmp79 + tmp178*tmp93 + tmp181*tmp93 + tmp195*tmp86 + tmp196*tmp203 + tmp196*tmp84 + tmp198*tmp199 + tmp202*tmp83 + tmp204*tmp79 + tmp207*tmp208 + tmp207*tmp90 + tmp208*tmp210 + tmp210*tmp90 + tmp211*tmp22 + 2*tmp212*tmp90 + tmp213*tmp214*tmp83 + tmp214*tmp215*tmp78 - tmp22*tmp37 - tmp76*(tmp46*tmp65 + tmp71*tmp74 + tmp75*(-tmp13*tmp38 + tmp21*tmp38)));
               aux_gfs[IDX4(MU1GF, i0, i1, i2)] = tmp230*(cf_dD0*(tmp158*tmp81 + tmp203*tmp228 + tmp208*tmp224 + tmp221*tmp222 + tmp222*tmp225 + tmp223*tmp83 + tmp224*tmp90 + tmp226*tmp227 + tmp85*tmp95) + tmp120*pow(tmp122, 2) + tmp122*tmp202 + tmp123*tmp220 + tmp141*tmp219 + tmp154*tmp98 + tmp157*tmp218 + tmp158*tmp175 + tmp158*tmp204 + tmp177*tmp84 + tmp181*tmp77*tmp83 + 2*tmp196*tmp226 + 2*tmp198*tmp219 + tmp200*tmp218 + tmp203*tmp213 + tmp207*tmp221 + tmp207*tmp225 + tmp208*tmp215 + 2*tmp210*tmp221 + tmp211*tmp83 + tmp212*tmp221 + tmp212*tmp225 + tmp213*tmp84 + tmp215*tmp90 + tmp229*tmp83 - tmp37*tmp83 - tmp76*(tmp217*tmp74 + tmp46*tmp75 + tmp65*(tmp121*tmp38 - tmp30*tmp38)));
               aux_gfs[IDX4(MU2GF, i0, i1, i2)] = tmp230*tmp9*(cf_dD0*(tmp158*tmp222 + tmp208*tmp228 + tmp221*tmp227 + tmp222*tmp232 + tmp223*tmp78 + tmp224*tmp233 + tmp224*tmp79 + tmp228*tmp90 + tmp231*tmp81) + pow(tmp117, 2)*tmp153 + tmp120*tmp219 + tmp122*tmp197*tmp201 + tmp124*tmp220 + tmp140*tmp218 + tmp157*tmp201*tmp98 + 2*tmp158*tmp207 + tmp158*tmp210 + tmp158*tmp212 + tmp175*tmp231 + tmp178*tmp89 + tmp180*tmp79 + tmp196*tmp221 + tmp196*tmp225 + tmp197*tmp218 + 2*tmp202*tmp98 + tmp208*tmp213 + tmp210*tmp232 + 2*tmp211*tmp78 + tmp212*tmp232 + tmp213*tmp90 + tmp215*tmp233 + tmp215*tmp79 + tmp229*tmp78 - tmp37*tmp78 - tmp76*(tmp217*tmp65 + tmp71*tmp75 + tmp74*(tmp116*tmp38 - tmp32*tmp38)));
            }
            
            
        } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++)
    } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++)
} // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++)
