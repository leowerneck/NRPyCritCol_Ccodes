// Part P1: Import needed header files
#include "BSSN_Ccodes/NGHOSTS.h" // A NRPy+-generated file, which is set based on FD_CENTDERIVS_ORDER.
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

// Part P2: Add needed #define's to set data type, the IDX4() macro, and the gridfunctions
// Part P2a: set REAL=double, so that all floating point numbers are stored to at least ~16 significant digits.
#define REAL double

// Step P3: Set free parameters for the numerical grid
const REAL RMAX = 64.0;

// Time coordinate parameters
REAL t_final;  /* Final time is set so that at t=t_final,
		* data at the origin have not been corrupted
		* by the approximate outer boundary condition */
const REAL CFL_FACTOR = 0.5; // Set the CFL Factor

// SinhSpherical parameters
REAL AMPL;
REAL SINHW;
const REAL REGRID_CRITERION = pow(200.0,1.0/30.0);
const int max_number_of_regrids = 30;

#define TFINAL (6.58)

// Step P3b: Free parameters for the spacetime evolution
REAL eta; // Gamma-driving shift condition parameter.
REAL diss_strength; // Kreiss-Oliger dissipation strenght parameter.

// Step P4: Implement the algorithm for upwinding.
//          *NOTE*: This upwinding is backwards from
//          usual upwinding algorithms, because the
//          upwinding control vector in BSSN (the shift)
//          acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) ( UpwindVecU > 0.0 ? 1.0 : 0.0 )

typedef struct __ID_inputs {
  int interp_stencil_size;
  int numlines_in_file;
  REAL *r_arr,*sf_arr,*psi4_arr,*alpha_arr;
} ID_inputs;

// Part P4b: Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//           data in a 1D array. In this case, consecutive values of "i"
//           (all other indices held to a fixed value) are consecutive in memory, where
//           consecutive values of "j" (fixing all other indices) are separated by
//           Nxx_plus_2NGHOSTS[0] elements in memory. Similarly, consecutive values of
//           "k" are separated by Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1] in memory, etc.
#define IDX4(g,i,j,k)							\
  ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * (k) ) )
// Assuming idx = IDX3(i,j,k). Much faster if idx can be reused over and over:
#define IDX4pt(g,idx)   ( (idx) + (Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2]) * (g) )

// Part P4c: Set #define's for BSSN gridfunctions. C code generated above
#include "BSSN_Ccodes/gridfunction_defines.h"

#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max)		\
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)

  // Step P5: Function for converting uniform grid coord
  //         (xx[0][i0],xx[1][i1],xx[2][i2]) to
  //          corresponding Cartesian coordinate.
void xxCart(REAL *xx[3],const int i0,const int i1,const int i2, REAL xCart[3]) {
  REAL xx0 = xx[0][i0];
  REAL xx1 = xx[1][i1];
  REAL xx2 = xx[2][i2];
#include "BSSN_Ccodes/xxCart.h"
}

// Step P6: Include basic functions needed to impose curvilinear
//          parity and boundary conditions.
#include "BSSN_Ccodes/curvilinear_parity_and_outer_boundary_conditions.h"

// Step P7: Function for enforcing the gammabar=gammahat constraint:
#include "BSSN_Ccodes/enforce_detgammabar_constraint.h"

// Step P8: Find the CFL-constrained timestep
REAL find_timestep(const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3],REAL *xx[3], const REAL CFL_FACTOR) {
  const REAL dxx0 = dxx[0], dxx1 = dxx[1], dxx2 = dxx[2];
  REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
  LOOP_REGION(NGHOSTS,Nxx_plus_2NGHOSTS[0]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[1]-NGHOSTS, NGHOSTS,Nxx_plus_2NGHOSTS[2]-NGHOSTS) {
    const REAL xx0 = xx[0][i0], xx1 = xx[1][i1];
    REAL ds_dirn0, ds_dirn1, ds_dirn2;
#include "BSSN_Ccodes/ds_dirn.h"
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
    // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
    dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
  }
  return dsmin*CFL_FACTOR;
}

// Part P9: Declare all functions for setting up SFC initial data.
/* Routines to interpolate the SFC solution and convert to ADM & T^{munu}: */
#include "ScalarField_Ccodes/scalarfield_interp.h"
#include "ScalarField_Ccodes/ID_scalar_field_ADM_quantities.h"
#include "ScalarField_Ccodes/ID_scalar_field_spherical.h"
#include "ScalarField_Ccodes/ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h"
#include "ScalarField_Ccodes/ID_scalarfield.h"

/* Next perform the basis conversion and compute all needed BSSN quantities */
#include "BSSN_Ccodes/ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs.h"
#include "BSSN_Ccodes/ID_BSSN__ALL_BUT_LAMBDAs.h"
#include "BSSN_Ccodes/ID_BSSN_lambdas.h"

// Part P10: Declare function for computing the Hamiltonian
//           constraint violation, which should converge to
//           zero with increasing numerical resolution.
void Hamiltonian_constraint(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3],
                            REAL *in_gfs, REAL *aux_gfs) {
#include "BSSN_Ccodes/Hamiltonian.h"
}

// Declare function for computing Momentum constraints violation.
void Momentum_constraints(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3],
                            REAL *in_gfs, REAL *aux_gfs) {
#include "BSSN_Ccodes/MomentumConstraint.h"
}

void compute_scalar_field_Tmunu(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3],
                            REAL *in_gfs, REAL *aux_gfs) {
#include "ScalarField_Ccodes/scalar_field_Tmunu.h"
}

// Part P11: Declare the function to evaluate the BSSN plus Scalar Field RHSs
void rhs_eval(const int Nxx[3],const int Nxx_plus_2NGHOSTS[3],const REAL dxx[3], REAL *xx[3],
              REAL *aux_gfs,const REAL *in_gfs,REAL *rhs_gfs) {
#include "BSSN_Ccodes/BSSN_plus_Scalar_Field_RHSs.h"
}

// Regridding functions
#include "Helpers_Ccodes/regrid.h"

// Output gridfunctions at the origin
#include "Helpers_Ccodes/output_central_values.h"

// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up scalar wave initial data
// Step 2: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
//         applying quadratic extrapolation outer boundary conditions.
// Step 3: Output relative error between numerical and exact solution.
// Step 4: Free all allocated memory
int main(int argc, const char *argv[]) {
    // Step 0a: Read command-line input, error out if nonconformant
    if(argc != 8 || atoi(argv[1]) < NGHOSTS) {
        fprintf(stderr,"Error: Expected three command-line arguments: ./SFC_Playground Nx0 Nx1 Nx2 sinhA sinhW diss_strength_eta diss_strength_KO,\n");
        fprintf(stderr,"where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
        fprintf(stderr,"Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
        exit(1);
    }
    // Step 0b: Set up numerical grid structure, first in space...
    if(atoi(argv[1])%2 != 0 || atoi(argv[2])%2 != 0 || atoi(argv[2])%2 != 0) {
        fprintf(stderr,"Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
        fprintf(stderr,"       For example, in case of angular directions, proper symmetry zones will not exist.\n");
        exit(1);
    }

    AMPL          = atof(argv[4]);
    SINHW         = atof(argv[5]);
    eta           = atof(argv[6]);
    diss_strength = atof(argv[7]);

    /* clear central values file */
    char central_values_filename_clear[100];
    sprintf(central_values_filename_clear,"out_central_values.dat");
    FILE *central_values_clear;
    central_values_clear = fopen(central_values_filename_clear,"w");
    fprintf(central_values_clear,"\n");
    fclose(central_values_clear);
    
    const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
    const int Nxx_plus_2NGHOSTS[3] = { Nxx[0]+2*NGHOSTS, Nxx[1]+2*NGHOSTS, Nxx[2]+2*NGHOSTS };
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1]*Nxx_plus_2NGHOSTS[2];
#include "BSSN_Ccodes/xxminmax.h"

    /* SFC INPUT ROUTINE */
    // Open the data file:
    char filename[100];
    sprintf(filename,"outputSFC.txt");
    FILE *in1Dsfc = fopen(filename, "r");
    if (in1Dsfc == NULL) {
      fprintf(stderr,"ERROR: could not open file %s\n",filename);
      exit(1);
    }
    // Count the number of lines in the data file:
    int numlines_in_file = count_num_lines_in_file(in1Dsfc);
    // Allocate space for all data arrays:
    REAL *r_arr     = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *sf_arr    = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *psi4_arr  = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    REAL *alpha_arr = (REAL *)malloc(sizeof(REAL)*numlines_in_file);
    
    // Read from the data file, filling in arrays
    if(read_datafile__set_arrays(in1Dsfc, r_arr,sf_arr,psi4_arr,alpha_arr) == 1) {
      fprintf(stderr,"ERROR WHEN READING FILE %s!\n",filename);
      exit(1);
    }
    fclose(in1Dsfc);

    ID_inputs SFC_in;

    const int interp_stencil_size = 2*NGHOSTS - 1;
    SFC_in.interp_stencil_size    = interp_stencil_size;
    SFC_in.numlines_in_file       = numlines_in_file;

    SFC_in.r_arr     = r_arr;
    SFC_in.sf_arr    = sf_arr;
    SFC_in.psi4_arr  = psi4_arr;
    SFC_in.alpha_arr = alpha_arr;
    /* END SFC INPUT ROUTINE */
    
    // Step 0c: Allocate memory for gridfunctions
    REAL *evol_gfs    = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    REAL *next_in_gfs = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    REAL *aux_gfs     = (REAL *)malloc(sizeof(REAL) * NUM_AUX_GFS  * Nxx_plus_2NGHOSTS_tot);
    REAL *k1_gfs      = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    REAL *k2_gfs      = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    REAL *k3_gfs      = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
    REAL *k4_gfs      = (REAL *)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);

    // Step 0d: Set up space and time coordinates
    // Step 0d.i: Set \Delta x^i on uniform grids.
    REAL dxx[3];
    for(int i=0;i<3;i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);

    // Step 0d.ii: Set up uniform coordinate grids
    REAL *xx[3];
    for(int i=0;i<3;i++) {
      xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
      for(int j=0;j<Nxx_plus_2NGHOSTS[i];j++) {
	xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxx[i]; // Cell-centered grid.
      }
    }
    
    // Step 0d.iii: Set timestep based on smallest proper distance between gridpoints and CFL factor
    REAL dt      = find_timestep(Nxx_plus_2NGHOSTS, dxx,xx, CFL_FACTOR);
    REAL t_final = MIN(TFINAL,AMPL);
    int N_final  = (int)(t_final / dt + 0.5); // The number of iterations in time. Add 0.5 to account for C rounding down integers.

    // Step 0e: Find ghostzone mappings and parities:
    gz_map *bc_gz_map = (gz_map *)malloc(sizeof(gz_map)*Nxx_plus_2NGHOSTS_tot);
    parity_condition *bc_parity_conditions = (parity_condition *)malloc(sizeof(parity_condition)*Nxx_plus_2NGHOSTS_tot);
    set_up_bc_gz_map_and_parity_conditions(Nxx_plus_2NGHOSTS,xx,dxx,xxmin,xxmax,  bc_gz_map, bc_parity_conditions);

    // Step 1: Convert the ADM initial data to BSSN initial data
    ID_BSSN__ALL_BUT_LAMBDAs(Nxx_plus_2NGHOSTS, xx, SFC_in, evol_gfs);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, NUM_EVOL_GFS, evol_gf_parity, evol_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
    ID_BSSN_lambdas(Nxx, Nxx_plus_2NGHOSTS, xx, dxx, evol_gfs);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, NUM_EVOL_GFS, evol_gf_parity, evol_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
    ID_scalarfield(Nxx_plus_2NGHOSTS, xx, SFC_in, evol_gfs);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions, NUM_EVOL_GFS, evol_gf_parity, evol_gfs);
    enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);

    compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);

    // Step 3a: Evaluate Hamiltonian constraint violation
    Hamiltonian_constraint(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
    
    // Step 3b: Evaluate Hamiltonian constraint violation
    Momentum_constraints(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
    apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);
    
    // Step 3: Start the timer, for keeping track of how fast the simulation is progressing.
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
    
    // Step 4: Integrate the initial data forward in time using the Method of Lines and RK4
    int n  = 0  , lapse_colapsed=0, regrid_count=0;
    REAL t = 0.0;
    
    while( t < t_final ) { // Main loop to progress forward in time.
                                 
      /***************************************************/
      /* Implement RK4 for Method of Lines timestepping: */
      /***************************************************/
      /* -= RK4: Step 1 of 4 =- */
      /* First evaluate k1 = RHSs expression             */
      rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, evol_gfs, k1_gfs);
      /* Next k1 -> k1*dt, and then set the input for    */
      /*    the next RHS eval call to y_n+k1/2           */
#pragma omp parallel for
      for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
	k1_gfs[i] *= dt;
	next_in_gfs[i] = evol_gfs[i] + k1_gfs[i]*0.5;
      }
      /* Finally, apply boundary conditions to           */
      /* next_in_gfs, so its data are set everywhere.    */
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, next_in_gfs);
      enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);
      compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, next_in_gfs, aux_gfs);
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);
                                 
      /* -= RK4: Step 2 of 4 =- */
      rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, next_in_gfs, k2_gfs);
#pragma omp parallel for
      for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
	k2_gfs[i] *= dt;
	next_in_gfs[i] = evol_gfs[i] + k2_gfs[i]*0.5;
      }
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, next_in_gfs);
      enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);
      compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, next_in_gfs, aux_gfs);
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);                       

      /* -= RK4: Step 3 of 4 =- */
      rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, next_in_gfs, k3_gfs);
#pragma omp parallel for
      for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
	k3_gfs[i] *= dt;
	next_in_gfs[i] = evol_gfs[i] + k3_gfs[i];
      }
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, next_in_gfs);
      enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, next_in_gfs);
      compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, next_in_gfs, aux_gfs);
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);
        
      /* -= RK4: Step 4 of 4 =- */
      rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, aux_gfs, next_in_gfs, k4_gfs);
#pragma omp parallel for
      for(int i=0;i<Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS;i++) {
	k4_gfs[i] *= dt;
	evol_gfs[i] = evol_gfs[i] + (1.0/6.0)*(k1_gfs[i] + 2.0*k2_gfs[i] + 2.0*k3_gfs[i] + k4_gfs[i]);
      }
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_EVOL_GFS,evol_gf_parity, evol_gfs);
      enforce_detgammabar_constraint(Nxx_plus_2NGHOSTS, xx, evol_gfs);
      compute_scalar_field_Tmunu(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, evol_gfs, aux_gfs);
      apply_bcs(Nxx, Nxx_plus_2NGHOSTS, bc_gz_map,bc_parity_conditions,NUM_AUX_GFS,aux_gf_parity, aux_gfs);

      // Check whether or not the lapse function has collapsed and output central values
      {
	const int i0 = NGHOSTS;
	const int i1 = NGHOSTS;
	const int i2 = NGHOSTS;
	const REAL alpha_i0p1  = evol_gfs[IDX4(ALPHAGF, i0,i1,i2)];
	const REAL alpha_i0p2  = evol_gfs[IDX4(ALPHAGF, i0+1,i1,i2)];
	const REAL alpha_i0p3  = evol_gfs[IDX4(ALPHAGF, i0+2,i1,i2)];
	const REAL alpha_c     = 3.0*alpha_i0p1 - 3.0*alpha_i0p2 + alpha_i0p3;

	const REAL sf_i0p1 = evol_gfs[IDX4(SFGF, i0,i1,i2)];
	const REAL sf_i0p2 = evol_gfs[IDX4(SFGF, i0+1,i1,i2)];
	const REAL sf_i0p3 = evol_gfs[IDX4(SFGF, i0+2,i1,i2)];
	const REAL sf_c    = 3.0*sf_i0p1    - 3.0*sf_i0p2    + sf_i0p3;

	char central_values_filename[100];
	sprintf(central_values_filename,"out_central_values.dat");
	FILE *central_values;
	central_values = fopen(central_values_filename,"a");
	fprintf(central_values,"%.15e %.15e %.15e\n",t,alpha_c,sf_c);
	fclose(central_values);

	if( alpha_c < 0.005 ) {
	  fprintf(stderr,"\nThe lapse function has collapsed...\n");
	  break;
	}
      }

#include "helpers_Ccodes/check_regrid_criterion.h"

      n++;
      t += dt;

      clock_gettime(CLOCK_REALTIME, &end);
      const long long unsigned int time_in_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
      const REAL s_per_iteration_avg          = ((REAL)time_in_ns / (REAL)n) / 1.0e9;
      const REAL execution_time_in_s          = (REAL)time_in_ns / 1.0e9;

      const REAL time_remaining        = t_final - t;
      const REAL time_per_hour         = (REAL)(dt * 3600.0 / s_per_iteration_avg);
      const REAL time_remaining_in_sec = (time_remaining / time_per_hour)*3600;

      const REAL num_RHS_pt_evals     = (REAL)(Nxx[0]*Nxx[1]*Nxx[2]) * 4.0 * (REAL)n; // 4 RHS evals per gridpoint for RK4
      const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / ((REAL)time_in_ns / 1.0e9);

      // Progress indicator printing to stderr
      fprintf(stderr,"%c[2K", 27); // Clear the line
      fprintf(stderr,"It: %d t=%.3f %.2f%% | Exec time: %.0f s; ETA %.0f s | t/h %.2f | gp/s %.2e | # of regrids: %d | sinhA = %.3lf | sinhW = %.3lf | eta_shift = %.2lf | KO = %.1lf\r",  // \r is carriage return, move cursor to the beginning of the line
      	      n, t, 100.0*t/t_final, execution_time_in_s,
      	      time_remaining_in_sec, time_per_hour, (REAL)RHS_pt_evals_per_sec, regrid_count,AMPL,SINHW,eta,diss_strength);
      fflush(stderr); // Flush the stderr buffer
        
    } // End main loop to progress forward in time.
    fprintf(stderr,"\n"); // Clear the line.

    /* Step 4: Free all allocated memory */
    free(bc_gz_map);
    free(bc_parity_conditions);
    free(aux_gfs);
    free(evol_gfs);
    free(k1_gfs);
    free(k2_gfs);
    free(k3_gfs);
    free(k4_gfs);
    for(int i=0;i<3;i++) free(xx[i]);
    return 0;

}
