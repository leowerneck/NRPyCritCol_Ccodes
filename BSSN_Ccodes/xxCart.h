/*
 *  Original SymPy expressions:
 *  "[xCart[0] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*cos(xx2)/(exp(1/SINHW) - exp(-1/SINHW)),
 *    xCart[1] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*sin(xx2)/(exp(1/SINHW) - exp(-1/SINHW)),
 *    xCart[2] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*cos(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
 */
{
   const double tmp0 = 1.0/SINHW;
   const double tmp1 = tmp0*xx0;
   const double tmp2 = AMPL*(exp(tmp1) - exp(-tmp1))/(exp(tmp0) - exp(-tmp0));
   const double tmp3 = tmp2*sin(xx1);
   xCart[0] = tmp3*cos(xx2);
   xCart[1] = tmp3*sin(xx2);
   xCart[2] = tmp2*cos(xx1);
}
