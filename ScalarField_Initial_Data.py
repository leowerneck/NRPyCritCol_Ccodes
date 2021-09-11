import numpy as np
import sys
from scipy.sparse import spdiags
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
    
def ScalarField_Initial_Data(outputname, P, R0, SIGMA, NR, RMAX, lapse_condition="Pre-collapsed"):
    ###################################################################
    # Part A.1: Setting up the initial condition for the scalar field #
    ###################################################################

    # Part A.1a: set up the necessary variables: number of radial points = NR
    #                                            the radial variable rr in [0,RMAX],
    #                                            \varphi(0,r) = ID_sf     (for scalar field),
    #                                            \Phi(0,r)    = ID_sf_dD0 (for the radial derivative of sf) ,
    #                                            \Pi(0,r)     = ID_sfM    (for scalar field conjutage Momentum).
    r  = np.linspace(0,RMAX,NR+1) # Set the r array
    dr = np.zeros(NR)
    for i in range(NR):
        dr[i] = r[1]-r[0]
    r  = np.delete(r-dr[0]/2,0)      # Shift the vector by -dr/2 and remove the negative entry

    dr2   = dr**2                    # Set the step size squared
    ID_sf = np.zeros(NR)

    ## Part A.1c: populating the varphi(0,r) array
    ID_sf = P * np.exp(-(r-R0)**2/SIGMA**2)

    ## Part A.1d: Create a function to evaluate Phi(r) := \partial_{r}\varphi
    def Phi(r):
        return -2 * P * (r - R0) * np.exp(-(r-R0)**2/SIGMA**2) / SIGMA**2

    # Set the main diagonal
    main_diag = np.pi * dr2 * Phi(r)**2 - 2

    # Update the first element of the main diagonal
    main_diag[0] += 1 - dr[0]/r[0]

    # Update the last element of the main diagonal
    main_diag[NR-1] += - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

    # Set the upper diagonal, ignoring the last point in the r array
    upper_diag = np.zeros(NR)
    upper_diag[1:] = 1 + dr[:-1]/r[:-1]

    # Set the lower diagonal, start counting the r array at the second element
    lower_diag = np.zeros(NR)
    lower_diag[:-1] = 1 - dr[1:]/r[1:]

    # Change the last term in the lower diagonal to its correct value
    lower_diag[NR-2] = 2

    # Set the sparse matrix A by adding up the three diagonals
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.spdiags.html
    A = spdiags([main_diag,upper_diag,lower_diag],[0,1,-1],NR,NR)

    # Then compress the sparse matrix A column wise, so that SciPy can invert it later
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    A = csc_matrix(A)

    # Set up the right-hand side of the linear system: s
    s = np.zeros(NR)

    # Update the last entry of the vector s
    s[NR-1] = - (2 * dr[NR-1] / r[NR-1])*(1 + dr[NR-1] / r[NR-1])

    # Compress the vector s column-wise
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csc_matrix.html
    s = csc_matrix(s)

    # Solve the sparse linear system using scipy
    # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.sparse.linalg.spsolve.html
    psi = spsolve(A, s.T)
    
    if lapse_condition == "Pre-collapsed":
        ID_alpha = psi**(-2)
        
        if sys.version_info[0] == 3:
            np.savetxt(outputname, list(zip( r, ID_sf, psi**4, ID_alpha )),fmt="%.15e")

        elif sys.version_info[0] == 2:
            np.savetxt(outputname, zip( r, ID_sf, psi**4, ID_alpha ),fmt="%.15e")

    elif lapse_condition == "Unity":
        ID_alpha = np.ones(NR)

        if sys.version_info[0] == 3:
            np.savetxt(outputname, list(zip( r, ID_sf, psi**4, ID_alpha )),fmt="%.15e")

        elif sys.version_info[0] == 2:
            np.savetxt(outputname, zip( r, ID_sf, psi**4, ID_alpha ),fmt="%.15e")
            
    else:
        print("Error: unknown lapse condition. Available options are: \"Pre-collapsed\" and \"Unity\"")
        return
        
    print("Generated the ADM initial data for the gravitational collapse \n" \
          "of a massless scalar field in Spherical coordinates.\n")
    print("Type of initial condition: Scalar field: \"Gaussian\" Shell\n"\
          "                         ADM quantities: Time-symmetric\n"\
          "                        Lapse condition: "+lapse_condition)
    print("Parameters: P       = "+str(P)+",\n" \
          "            R0      = "+str(R0)+",\n"   \
          "            SIGMA   = "+str(SIGMA)+",\n"   \
          "            RMAX    = "+str(RMAX)+",\n"   \
          "            NR      = "+str(NR)+",\n"
          "            Outfile = "+str(outputname)+".\n")

if __name__ == '__main__':

    if len(sys.argv) != 7:
        print("ERROR: correct usage is:")
        print("python ScalarField_Initial_Data.py <output_filename> <pulse_amplitude> <pulse_center> <pulse_width> <Nr> <R_max>")
        sys.exit(1)

    output_filename = str(  sys.argv[1])
    pulse_amplitude = float(sys.argv[2])
    pulse_center    = float(sys.argv[3])
    pulse_width     = float(sys.argv[4])
    Nr              = int(  sys.argv[5])
    R_max           = float(sys.argv[6])
    
    ScalarField_Initial_Data(output_filename,
                             pulse_amplitude,
                             pulse_center,
                             pulse_width,
                             Nr,
                             R_max)
