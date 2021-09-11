const REAL invdx0 = 1.0/dxx[0];
#pragma omp parallel for
for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++) {
    for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++) {
        const REAL xx1 = xx[1][i1];
        for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++) {
            const REAL xx0 = xx[0][i0];
            {
               /* 
                * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                */
               const double hDD00 = in_gfs[IDX4(HDD00GF, i0,i1,i2)];
               const double hDD01 = in_gfs[IDX4(HDD01GF, i0,i1,i2)];
               const double hDD02 = in_gfs[IDX4(HDD02GF, i0,i1,i2)];
               const double hDD11 = in_gfs[IDX4(HDD11GF, i0,i1,i2)];
               const double hDD12 = in_gfs[IDX4(HDD12GF, i0,i1,i2)];
               const double hDD22 = in_gfs[IDX4(HDD22GF, i0,i1,i2)];
               const double vetU0 = in_gfs[IDX4(VETU0GF, i0,i1,i2)];
               const double vetU1 = in_gfs[IDX4(VETU1GF, i0,i1,i2)];
               const double vetU2 = in_gfs[IDX4(VETU2GF, i0,i1,i2)];
               const double cf = in_gfs[IDX4(CFGF, i0,i1,i2)];
               const double alpha = in_gfs[IDX4(ALPHAGF, i0,i1,i2)];
               const double sf_i0m2_i1_i2 = in_gfs[IDX4(SFGF, i0-2,i1,i2)];
               const double sf_i0m1_i1_i2 = in_gfs[IDX4(SFGF, i0-1,i1,i2)];
               const double sf_i0p1_i1_i2 = in_gfs[IDX4(SFGF, i0+1,i1,i2)];
               const double sf_i0p2_i1_i2 = in_gfs[IDX4(SFGF, i0+2,i1,i2)];
               const double sfM = in_gfs[IDX4(SFMGF, i0,i1,i2)];
               const double sf_dD0 = invdx0*(-2.0/3.0*sf_i0m1_i1_i2 + (1.0/12.0)*sf_i0m2_i1_i2 + (2.0/3.0)*sf_i0p1_i1_i2 - 1.0/12.0*sf_i0p2_i1_i2);
               /* 
                * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                */
               const double tmp0 = pow(alpha, -2);
               const double tmp1 = pow(sfM, 2);
               const double tmp2 = exp(8*cf);
               const double tmp3 = 1.0/SINHW;
               const double tmp4 = exp(tmp3) - exp(-tmp3);
               const double tmp5 = pow(AMPL, 4)/pow(tmp4, 4);
               const double tmp6 = tmp2*tmp5;
               const double tmp7 = sin(xx1);
               const double tmp8 = pow(tmp7, 2);
               const double tmp9 = tmp3*xx0;
               const double tmp10 = exp(tmp9);
               const double tmp11 = exp(-tmp9);
               const double tmp12 = tmp10 - tmp11;
               const double tmp13 = pow(tmp12, 4)*tmp8;
               const double tmp14 = pow(hDD12, 2)*tmp13;
               const double tmp15 = pow(tmp12, 2);
               const double tmp16 = pow(AMPL, 2);
               const double tmp17 = pow(tmp4, 2);
               const double tmp18 = tmp16/tmp17;
               const double tmp19 = tmp15*tmp18;
               const double tmp20 = hDD11*tmp19 + tmp19;
               const double tmp21 = tmp19*tmp8;
               const double tmp22 = hDD22*tmp21 + tmp21;
               const double tmp23 = tmp2*tmp22;
               const double tmp24 = tmp10*tmp3 + tmp11*tmp3;
               const double tmp25 = pow(tmp24, 2);
               const double tmp26 = exp(12*cf);
               const double tmp27 = hDD01*hDD02;
               const double tmp28 = tmp18*tmp25;
               const double tmp29 = hDD00*tmp28 + tmp28;
               const double tmp30 = tmp26*tmp5;
               const double tmp31 = tmp15*tmp25;
               const double tmp32 = pow(hDD02, 2)*tmp31;
               const double tmp33 = tmp22*tmp26;
               const double tmp34 = pow(hDD01, 2)*tmp31;
               const double tmp35 = tmp20*tmp29;
               const double tmp36 = 1.0/(2*pow(AMPL, 6)*hDD12*tmp13*tmp25*tmp26*tmp27/pow(tmp4, 6) - tmp14*tmp29*tmp30 - tmp20*tmp30*tmp32*tmp8 - tmp33*tmp34*tmp5 + tmp33*tmp35);
               const double tmp37 = tmp36*(-tmp14*tmp6 + tmp20*tmp23);
               const double tmp38 = (1.0/2.0)*pow(sf_dD0, 2)*tmp37 - 1.0/2.0*tmp1;
               const double tmp39 = tmp0*tmp38;
               const double tmp40 = vetU0/tmp24;
               const double tmp41 = sfM/alpha;
               const double tmp42 = tmp4/AMPL;
               const double tmp43 = tmp41*tmp42;
               const double tmp44 = sf_dD0*tmp37 - tmp40*tmp43;
               const double tmp45 = tmp39*tmp42;
               const double tmp46 = 1.0/tmp12;
               const double tmp47 = tmp46*vetU1;
               const double tmp48 = pow(tmp12, 3);
               const double tmp49 = tmp6*tmp8;
               const double tmp50 = hDD02*tmp24;
               const double tmp51 = hDD01*tmp24;
               const double tmp52 = tmp12*tmp18;
               const double tmp53 = tmp36*(hDD12*tmp48*tmp49*tmp50 - tmp23*tmp51*tmp52);
               const double tmp54 = sf_dD0*tmp53 - tmp43*tmp47;
               const double tmp55 = vetU2/tmp7;
               const double tmp56 = tmp46*tmp55;
               const double tmp57 = hDD12*tmp7;
               const double tmp58 = tmp36*(-tmp2*tmp20*tmp50*tmp52*tmp7 + tmp48*tmp51*tmp57*tmp6);
               const double tmp59 = sf_dD0*tmp58 - tmp43*tmp56;
               const double tmp60 = tmp0*tmp17/tmp16;
               const double tmp61 = tmp40*tmp60;
               const double tmp62 = tmp60/tmp15;
               aux_gfs[IDX4(T4UU00GF, i0, i1, i2)] = tmp0*tmp1 + tmp39;
               aux_gfs[IDX4(T4UU01GF, i0, i1, i2)] = -tmp40*tmp45 + tmp41*tmp44;
               aux_gfs[IDX4(T4UU02GF, i0, i1, i2)] = tmp41*tmp54 - tmp45*tmp47;
               aux_gfs[IDX4(T4UU03GF, i0, i1, i2)] = tmp41*tmp59 - tmp45*tmp56;
               aux_gfs[IDX4(T4UU11GF, i0, i1, i2)] = -tmp38*(tmp37 - tmp60*pow(vetU0, 2)/tmp25) + pow(tmp44, 2);
               aux_gfs[IDX4(T4UU12GF, i0, i1, i2)] = -tmp38*(-tmp47*tmp61 + tmp53) + tmp44*tmp54;
               aux_gfs[IDX4(T4UU13GF, i0, i1, i2)] = -tmp38*(-tmp56*tmp61 + tmp58) + tmp44*tmp59;
               aux_gfs[IDX4(T4UU22GF, i0, i1, i2)] = -tmp38*(tmp36*(tmp23*tmp29 - tmp32*tmp49) - tmp62*pow(vetU1, 2)) + pow(tmp54, 2);
               aux_gfs[IDX4(T4UU23GF, i0, i1, i2)] = -tmp38*(tmp36*(-tmp19*tmp2*tmp29*tmp57 + tmp27*tmp31*tmp6*tmp7) - tmp55*tmp62*vetU1) + tmp54*tmp59;
               aux_gfs[IDX4(T4UU33GF, i0, i1, i2)] = -tmp38*(tmp36*(tmp2*tmp35 - tmp34*tmp6) - tmp62*pow(vetU2, 2)/tmp8) + pow(tmp59, 2);
            }
            
            
        } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++)
    } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++)
} // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++)
