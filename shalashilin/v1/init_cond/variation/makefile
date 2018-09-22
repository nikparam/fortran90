F = f2py

FLAGS = -c -m

SOURCE = ../source/hmatrix_legendre.f90 \
	 ../../D_coeff/source/smatrix.f90 \
	 ../../D_coeff/source/change_var_mod.f90 \
	 ../../../../legendre_polynomials/legendre_polynomial.f90 \
	 ../../potential/pot_der.f90

NAME = fort_subroutines

all: clean

	@$(F) $(FLAGS) $(NAME) $(SOURCE)
	@mv $(NAME)* $(NAME).so

clean: 

	rm -rf *.so *.pyf
