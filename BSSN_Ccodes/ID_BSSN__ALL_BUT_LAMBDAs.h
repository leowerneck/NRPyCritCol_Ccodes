void ID_BSSN__ALL_BUT_LAMBDAs(const int Nxx_plus_2NGHOSTS[3],REAL *xx[3],ID_inputs other_inputs,REAL *in_gfs) {
#pragma omp parallel for
for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++) {
    const REAL xx2 = xx[2][i2];
    for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++) {
        const REAL xx1 = xx[1][i1];
        for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++) {
            const REAL xx0 = xx[0][i0];
            const int idx = IDX3(i0,i1,i2);
            const REAL xx0xx1xx2[3] = {xx0,xx1,xx2};
            ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs(xx0xx1xx2,other_inputs,
                                &in_gfs[IDX4pt(HDD00GF,idx)],&in_gfs[IDX4pt(HDD01GF,idx)],&in_gfs[IDX4pt(HDD02GF,idx)],
                                &in_gfs[IDX4pt(HDD11GF,idx)],&in_gfs[IDX4pt(HDD12GF,idx)],&in_gfs[IDX4pt(HDD22GF,idx)],
                                &in_gfs[IDX4pt(ADD00GF,idx)],&in_gfs[IDX4pt(ADD01GF,idx)],&in_gfs[IDX4pt(ADD02GF,idx)],
                                &in_gfs[IDX4pt(ADD11GF,idx)],&in_gfs[IDX4pt(ADD12GF,idx)],&in_gfs[IDX4pt(ADD22GF,idx)],
                                &in_gfs[IDX4pt(TRKGF,idx)],
                                &in_gfs[IDX4pt(VETU0GF,idx)],&in_gfs[IDX4pt(VETU1GF,idx)],&in_gfs[IDX4pt(VETU2GF,idx)],
                                &in_gfs[IDX4pt(BETU0GF,idx)],&in_gfs[IDX4pt(BETU1GF,idx)],&in_gfs[IDX4pt(BETU2GF,idx)],
                                &in_gfs[IDX4pt(ALPHAGF,idx)],&in_gfs[IDX4pt(CFGF,idx)]);
            
        } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++)
    } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++)
} // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++)
}
