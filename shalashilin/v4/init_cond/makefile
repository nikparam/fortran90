F = gfortran 

DIR = ./old
FILES = *.o *.mod *.out *.tmp
EXEC = ./main

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

main_get_grid: ./source/add_to_list.f90  ./source/get_grid.f90 ./source/get_states.f90 ./source/main.f90  ../potential/pot_der.f90 ../D_coeff/source/smatrix.f90 ../D_coeff/source/hmatrix.f90 ../D_coeff/source/change_var_mod.f90 ../D_coeff/source/normalization.f90

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC)

main_get_hollow_grid: ./source/add_to_list.f90  ./source/get_hollow_grid.f90 ./source/get_states.f90 ./source/main.f90  ../potential/pot_der.f90 ../D_coeff/source/smatrix.f90 ../D_coeff/source/hmatrix.f90 ../D_coeff/source/change_var_mod.f90 ../D_coeff/source/normalization.f90

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC)

main_get_cross_grid: ./source/add_to_list.f90  ./source/get_cross_grid.f90 ./source/get_states.f90 ./source/main.f90  ../potential/pot_der.f90 ../D_coeff/source/smatrix.f90 ../D_coeff/source/hmatrix.f90 ../D_coeff/source/change_var_mod.f90 ../D_coeff/source/normalization.f90

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC)

clean:

	@rm -rf $(FILES)
