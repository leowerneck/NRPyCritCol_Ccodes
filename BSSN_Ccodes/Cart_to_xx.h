/*
 *  Original SymPy expressions:
 *  "[Cart_to_xx0_inbounds = SINHW*asinh(sqrt(Cartx**2 + Carty**2 + Cartz**2)*sinh(1/SINHW)/AMPL),
 *    Cart_to_xx1_inbounds = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2)),
 *    Cart_to_xx2_inbounds = atan2(Carty, Cartx)]"
 */
{
   const double tmp0 = sqrt(pow(Cartx, 2) + pow(Carty, 2) + pow(Cartz, 2));
   Cart_to_xx0_inbounds = SINHW*asinh(tmp0*sinh(1.0/SINHW)/AMPL);
   Cart_to_xx1_inbounds = acos(Cartz/tmp0);
   Cart_to_xx2_inbounds = atan2(Carty, Cartx);
}
