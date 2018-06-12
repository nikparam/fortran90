F = mpif90

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 
EXEC = main_pot	
SOURCE = ./source/main.f90 \
	./source/trajectories.f90 \
	./source/externalqp.f90 \
	../../../odepack/vode.f \
	../../../odepack/dgesl.f \
	../../../odepack/dgefa.f \
	../../../odepack/dgbfa.f \
	../../../odepack/dgbsl.f

main_pot: $(SOURCE) ../potential/pot_der.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_morse: $(SOURCE) ../potential/morse_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_exp: $(SOURCE) ../potential/exp_potential.f90

	@$(F) $^ $(LD_FLAGS) -o $(EXEC)

clean:

	@rm -rf *.out *.o $(EXEC)
