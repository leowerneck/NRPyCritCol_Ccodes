/*
 *  Original SymPy expressions:
 *  "[ds_dirn0 = AMPL*dxx0*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)/(exp(1/SINHW) - exp(-1/SINHW)),
 *    ds_dirn1 = AMPL*dxx1*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1/SINHW) - exp(-1/SINHW)),
 *    ds_dirn2 = AMPL*dxx2*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
 */
{
   const double tmp0 = 1.0/SINHW;
   const double tmp1 = tmp0*xx0;
   const double tmp2 = exp(tmp1);
   const double tmp3 = exp(-tmp1);
   const double tmp4 = AMPL/(exp(tmp0) - exp(-tmp0));
   const double tmp5 = tmp4*(tmp2 - tmp3);
   ds_dirn0 = dxx0*tmp4*(tmp0*tmp2 + tmp0*tmp3);
   ds_dirn1 = dxx1*tmp5;
   ds_dirn2 = dxx2*tmp5*sin(xx1);
}
