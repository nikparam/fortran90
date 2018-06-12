F = gfortran

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

EXEC = ./main

SOURCE =  ./source/external_D_mod.f90 \
          ./source/coeff_matrix.f90  \
          ./source/change_var_mod.f90 \
          ./source/smatrix.f90 \
          ./source/sinvert.f90 \
          ./source/norm_energy.f90 \
          ./source/normalization.f90 \
          ../../../odepack/zvode.f \
          ../../../odepack/zgesl.f \
          ../../../odepack/zgefa.f \
          ../../../odepack/zgbfa.f \
          ../../../odepack/zgbsl.f 

legendre = ./source/main_num_l.f90 \
           ./source/hmatrix_legendre.f90 \
	   ../../../legendre_polynomials/legendre_polynomial.f90 

hermite = ./source/main_num_h.f90 \
	  ./source/hmatrix_hermite.f90 \
	  ../../../hermite_polynomials/hermite_polynomial.f90 

analytical = ./source/main_analytical.f90 \
	     ./source/hmatrix.f90


main_pot_num_l:  $(SOURCE) $(legendre) ../potential/pot_der.f90 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_pot_num_h: $(SOURCE) $(hermite) ../potential/pot_der.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_pot_an: $(SOURCE) $(analytical) ../potential/pot_der.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_morse_num: $(SOURCE) $(legendre) ../potential/morse_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_exp_num: $(SOURCE) $(legendre) ../potential/exp_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

clean:

	@rm -rf *.out *.o $(EXEC)

