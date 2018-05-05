F = gfortran 

DIR = ./old
DIR_sym = ./old/sym
DIR_svd = ./old/svd
DIR_chl = ./old/chl
FILES = *.o *.mod *.out *.tmp
EXEC = ./main*
EXEC_sym = ./main_lowdin
EXEC_svd = ./main_svd
EXEC_chl = ./main_chl

SOURCE = ./source/add_to_list.f90 \
         ./source/get_cross_grid.f90 \
         ./source/norm_energy.f90 \
         ./source/main.f90  \
         ../potential/morse_potential.f90 \
         ../D_coeff/source/smatrix.f90 \
         ../D_coeff/source/change_var_mod.f90 \
         ../D_coeff/source/normalization.f90 


legendre = ./source/hmatrix_legendre.f90 \
 	   ../../../legendre_polynomials/legendre_polynomial.f90

lowdin = ./source/get_states_sym.f90 
svd = ./source/get_states_svd.f90 
chl = ./source/get_states_chl.f90 

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 

all: clean cg_low cg_svd cg_chl

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

cg_low: $(SOURCE) $(legendre) $(lowdin)

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_sym)

cg_svd: $(SOURCE) $(legendre) $(svd)	

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_svd)

cg_chl:	$(SOURCE) $(legendre) $(chl)	

	@$(F)$ $^ $(LD_FLAGS) -o $(EXEC_chl)

clean:

	@rm -rf $(FILES) $(DIR) $(EXEC) 
