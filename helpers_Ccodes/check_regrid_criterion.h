  {
    if( ( t > 5.100 ) && ( regrid_count == 0 ) ) {

	diss_strength = 1.0;
  
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.19, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);
      	regrid_count++;

      }

      if( ( t > 5.125 ) && ( regrid_count == 1 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.18, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.150 ) && ( regrid_count == 2 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.17, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.175 ) && ( regrid_count == 3 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.16, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.200 ) && ( regrid_count == 4 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.15, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.225 ) && ( regrid_count == 5 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.14, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.250 ) && ( regrid_count == 6 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.13, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.275 ) && ( regrid_count == 7 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.12, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.300 ) && ( regrid_count == 8 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.11, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 5.325 ) && ( regrid_count == 9 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_SINHW( regrid_stencil_size, 0.10, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.46 ) && ( regrid_count == 10 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,48, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.47 ) && ( regrid_count == 11 ) ) {

      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,32, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.48 ) && ( regrid_count == 12 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,24, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.49 ) && ( regrid_count == 13 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,16, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.552 ) && ( regrid_count == 14 ) ) {

        diss_strength = 7.0;
        
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
	regrid_SINHW( regrid_stencil_size,0.09, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.553 ) && ( regrid_count == 15 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,14, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.554 ) && ( regrid_count == 16 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,12, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.555 ) && ( regrid_count == 17 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,10, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }

      if( ( t > 6.556 ) && ( regrid_count == 18 ) ) {
	
      	const int regrid_stencil_size = 2*NGHOSTS - 1;
      	regrid_AMPL( regrid_stencil_size,8, Nxx_plus_2NGHOSTS, xx, next_in_gfs, evol_gfs );
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      	enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      	compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      	apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      	dt = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);

      	// Adjust t_final
      	t_final = MIN(t_final,t+AMPL);

      	N_final = (int)( t_final / dt + 0.5);

      	regrid_count++;

      }
  }
