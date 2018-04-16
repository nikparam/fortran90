F = gfortran

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

EXEC = ./main
	    
main_pot_num_l: ./source/main_num_l.f90 \
		 ./source/external_D_mod.f90 \
		 ./source/coeff_matrix.f90  \
		 ./source/change_var_mod.f90 \
		 ./source/smatrix.f90 \
		 ./source/sinvert.f90 \
		 ./source/norm_energy.f90 \
		 ./source/normalization.f90 \
		 ../potential/pot_der.f90 \
		 ../../../odepack/zvode.f \
		 ../../../odepack/zgesl.f \
		 ../../../odepack/zgefa.f \
		 ../../../odepack/zgbfa.f \
		 ../../../odepack/zgbsl.f \
		 ./source/hmatrix_legendre.f90 \
		 ../../../legendre_polynomials/legendre_polynomial.f90 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_pot_num_h: ./source/main_num_h.f90 \
		 ./source/external_D_mod.f90 \
		 ./source/coeff_matrix.f90  \
		 ./source/change_var_mod.f90 \
		 ./source/smatrix.f90 \
		 ./source/sinvert.f90 \
		 ./source/norm_energy.f90 \
		 ./source/normalization.f90 \
		 ../potential/pot_der.f90 \
		 ../../../odepack/zvode.f \
		 ../../../odepack/zgesl.f \
		 ../../../odepack/zgefa.f \
		 ../../../odepack/zgbfa.f \
		 ../../../odepack/zgbsl.f \
		 ./source/hmatrix_hermite.f90 \
		 ../../../hermite_polynomials/hermite_polynomial.f90 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_pot_an: ./source/main_analytical.f90 \
		 ./source/external_D_mod.f90 \
		 ./source/coeff_matrix.f90  \
		 ./source/change_var_mod.f90 \
		 ./source/smatrix.f90 \
		 ./source/hmatrix.f90 \
		 ./source/sinvert.f90 \
		 ./source/norm_energy.f90 \
		 ./source/normalization.f90 \
		 ../potential/pot_der.f90 \
		 ../../../odepack/zvode.f \
		 ../../../odepack/zgesl.f \
		 ../../../odepack/zgefa.f \
		 ../../../odepack/zgbfa.f \
		 ../../../odepack/zgbsl.f 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)


main_pot_morse: ./source/main.f90 ./source/external_D.f90 ./source/change_var.f90 ./source/smatrix.f90 ./source/hmatrix.f90 ./source/sinvert.f90 ./source/norm_energy.f90 ./source/normalization.f90 ../potential/morse_potential.f90 ../../../odepack/zvode.f ../../../odepack/zgesl.f ../../../odepack/zgefa.f ../../../odepack/zgbfa.f ../../../odepack/zgbsl.f 

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

clean:

	@rm -rf *.out *.o $(EXEC)

