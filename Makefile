CC     = gcc
CFLAGS = -fopenmp -march=native -Ofast
CLIBS  = -lm

SRC    = NRPyCritCol.c

BSSN_INCLUDE_FILES = BSSN_plus_Scalar_Field_RHSs.h                         \
                     Cart_to_xx.h                                          \
                     Hamiltonian.h                                         \
                     ID_ADM_xx0xx1xx2_to_BSSN_xx0xx1xx2__ALL_BUT_LAMBDAs.h \
                     ID_BSSN__ALL_BUT_LAMBDAs.h                            \
                     ID_BSSN_lambdas.h                                     \
                     MomentumConstraint.h                                  \
                     NGHOSTS.h                                             \
                     curvilinear_parity_and_outer_boundary_conditions.h    \
                     ds_dirn.h                                             \
                     enforce_detgammabar_constraint.h                      \
                     gridfunction_defines.h                                \
                     set_parity_conditions.h                               \
                     xxCart.h                                              \
                     xxminmax.h

SCALARFIELD_INCLUDE_FILES = ID_scalar_field_ADM_quantities.h               \
                            ID_scalar_field_spherical.h                    \
                            ID_scalarfield.h                               \
                            ID_scalarfield_xx0xx1xx2_to_BSSN_xx0xx1xx2.h   \
                            scalar_field_Tmunu.h                           \
                            scalarfield_interp.h

HELPERS_INCLUDE_FILES = check_regrid_criterion.h                           \
                        output_central_values.h                            \
                        regrid.h


OBJ = $(SRC:.c=.o)
EXE = $(SRC:.c=)
INC = $(addprefix BSSN_Ccodes/,       $(BSSN_INCLUDE_FILES)       ) \
      $(addprefix ScalarField_Ccodes/,$(SCALARFIELD_INCLUDE_FILES)) \
      $(addprefix helpers_Ccodes/,    $(HELPERS_INCLUDE_FILES)    )

all: $(EXE) runscript.sh

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $< -o $@ $(CLIBS)

$(OBJ): %.o : %.c $(INC)
	$(CC) $(CFLAGS) -c $<

runscript.sh: generate_runscript.py
	python generate_runscript.py

clean:
	rm -f $(OBJ) $(EXE) runscript.sh

veryclean: clean
	rm -f *.dat *.txt *.png
