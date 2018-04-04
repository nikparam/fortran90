F = gfortran 

DIR = ./old
DIR_sym = ./old/sym
DIR_svd = ./old/svd
DIR_chl = ./old/chl
FILES = *.o *.mod *.out *.tmp
EXEC = ./main*
EXEC_sym = ./main_sym
EXEC_svd = ./main_svd
EXEC_chl = ./main_chl

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

all: clean cg_sym cg_svd cg_chl

	@mkdir $(DIR)
	@$(EXEC_sym) 
	@mkdir $(DIR_sym)
	@mv *.out $(DIR_sym)

	@$(EXEC_svd) 
	@mkdir $(DIR_svd)
	@mv *.out $(DIR_svd)

	@$(EXEC_chl) 
	@mkdir $(DIR_chl)
	@mv *.out $(DIR_chl)

cg_sym: ./source/add_to_list.f90 \
        ./source/get_cross_grid.f90 \
        ./source/get_states_sym.f90 \
        ./source/main.f90  \
        ../potential/pot_der.f90 \
        ../D_coeff/source/smatrix.f90 \
        ../D_coeff/source/hmatrix.f90 \
        ../D_coeff/source/change_var_mod.f90 \
        ../D_coeff/source/normalization.f90 \
        ../D_coeff/source/norm_energy.f90 \

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_sym)

cg_svd: ./source/add_to_list.f90 \
        ./source/get_cross_grid.f90 \
        ./source/get_states_svd.f90 \
        ./source/main.f90  \
        ../potential/pot_der.f90 \
        ../D_coeff/source/smatrix.f90 \
        ../D_coeff/source/hmatrix.f90 \
        ../D_coeff/source/change_var_mod.f90 \
        ../D_coeff/source/normalization.f90 \
        ../D_coeff/source/norm_energy.f90 \
        ../D_coeff/source/sinvert.f90

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_svd)

cg_chl: ./source/add_to_list.f90 \
        ./source/get_cross_grid.f90 \
        ./source/get_states_chl.f90 \
        ./source/main.f90  \
        ../potential/pot_der.f90 \
        ../D_coeff/source/smatrix.f90 \
        ../D_coeff/source/hmatrix.f90 \
        ../D_coeff/source/change_var_mod.f90 \
        ../D_coeff/source/normalization.f90 \
        ../D_coeff/source/norm_energy.f90 \

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_chl)

clean:

	@rm -rf $(FILES) $(DIR) $(EXEC) 
