F = gfortran

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

EXEC = ./main

SOURCE =  ./source_new/main.f90 \
	  ./source_new/external_D.f90 \
          ./source_new/coeff_matrix.f90  \
          ./source_new/change_var.f90 \
          ./source_new/smatrix.f90 \
	  ./source_new/hmatrix.f90 \
          ./source_new/sinvert.f90 \
          ./source_new/normalization.f90 \
          ../../../odepack/zvode.f \
          ../../../odepack/zgesl.f \
          ../../../odepack/zgefa.f \
          ../../../odepack/zgbfa.f \
          ../../../odepack/zgbsl.f \
	  ../../../legendre_polynomials/legendre_polynomial.f90 
#	  ../../../../hermite_polynomials/hermite_polynomial.f90 

main_pot:  $(SOURCE) ../potential/pot_der.f90 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_morse: $(SOURCE) ../potential/morse_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_exp: $(SOURCE) ../potential/exp_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

clean:

	@rm -rf *.out *.o $(EXEC)

