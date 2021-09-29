void enforce_detgammabar_constraint(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3], REAL *in_gfs) {

for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++) {
    for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++) {
        const REAL xx1 = xx[1][i1];
#pragma omp parallel for
        for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++) {
            const REAL xx0 = xx[0][i0];
            /* 
             * NRPy+ Finite Difference Code Generation, Step 1 of 1: Read from main memory and compute finite difference stencils:
             */
            const double hDD00 = in_gfs[IDX4(HDD00GF, i0,i1,i2)];
            const double hDD01 = in_gfs[IDX4(HDD01GF, i0,i1,i2)];
            const double hDD02 = in_gfs[IDX4(HDD02GF, i0,i1,i2)];
            const double hDD11 = in_gfs[IDX4(HDD11GF, i0,i1,i2)];
            const double hDD12 = in_gfs[IDX4(HDD12GF, i0,i1,i2)];
            const double hDD22 = in_gfs[IDX4(HDD22GF, i0,i1,i2)];
            /* 
             * NRPy+ Finite Difference Code Generation, Step 2 of 1: Evaluate SymPy expressions and write to main memory:
             */
            const double tmp0 = pow(AMPL, 2);
            const double tmp1 = sin(xx1);
            const double tmp2 = pow(tmp1, 2);
            const double tmp3 = 1.0/SINHW;
            const double tmp4 = exp(tmp3) - exp(-tmp3);
            const double tmp5 = tmp3*xx0;
            const double tmp6 = exp(tmp5);
            const double tmp7 = exp(-tmp5);
            const double tmp8 = tmp6 - tmp7;
            const double tmp9 = pow(tmp8, 4);
            const double tmp10 = pow(tmp3*tmp6 + tmp3*tmp7, 2);
            const double tmp11 = tmp10*tmp9/pow(tmp4, 6);
            const double tmp12 = tmp0/pow(tmp4, 2);
            const double tmp13 = tmp10*tmp12;
            const double tmp14 = hDD00*tmp13 + tmp13;
            const double tmp15 = pow(AMPL, 4)/pow(tmp4, 4);
            const double tmp16 = pow(tmp8, 2);
            const double tmp17 = tmp12*tmp16;
            const double tmp18 = hDD11*tmp17 + tmp17;
            const double tmp19 = tmp10*tmp15*tmp16;
            const double tmp20 = tmp17*tmp2;
            const double tmp21 = hDD22*tmp20 + tmp20;
            const double tmp22 = tmp0*cbrt(1.0/(2*pow(AMPL, 6)*hDD01*hDD02*hDD12*tmp11*tmp2 - pow(hDD01, 2)*tmp19*tmp21 - pow(hDD02, 2)*tmp18*tmp19*tmp2 - pow(hDD12, 2)*tmp14*tmp15*tmp2*tmp9 + tmp14*tmp18*tmp21))*pow(fabs(tmp1), 2.0/3.0)*cbrt(fabs(tmp11));
            in_gfs[IDX4(HDD00GF, i0, i1, i2)] = tmp22*(hDD00 + 1) - 1;
            in_gfs[IDX4(HDD01GF, i0, i1, i2)] = hDD01*tmp22;
            in_gfs[IDX4(HDD02GF, i0, i1, i2)] = hDD02*tmp22;
            in_gfs[IDX4(HDD11GF, i0, i1, i2)] = tmp22*(hDD11 + 1) - 1;
            in_gfs[IDX4(HDD12GF, i0, i1, i2)] = hDD12*tmp22;
            in_gfs[IDX4(HDD22GF, i0, i1, i2)] = tmp22*(hDD22 + 1) - 1;
            
            
        } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++)
    } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++)
} // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++)
}
