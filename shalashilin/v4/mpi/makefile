F = mpif90

LD_FLAGS = -L/usr/lib/ -lblas -L/usr/lib/ -llapack 
EXEC = main_pot	
				    
main_pot: ./source/main.f90 ./source/trajectories.f90 ./source/externalqp.f90 \
	      ../potential/pot_der.f90 \
	      ../../../odepack/vode.f ../../../odepack/dgesl.f ../../../odepack/dgefa.f \
				      ../../../odepack/dgbfa.f ../../../odepack/dgbsl.f

	$(F) $^ $(LD_FLAGS) -o $(EXEC)

main_pot_morse: ./source/main.f90 ./source/trajectories.f90 ./source/externalqp.f90 \
	        ../potential/morse_potential.f90 \
	        ../../../odepack/vode.f ../../../odepack/dgesl.f ../../../odepack/dgefa.f \
				      ../../../odepack/dgbfa.f ../../../odepack/dgbsl.f

	$(F) $^ $(LD_FLAGS) -o $(EXEC)

clean: 

	@rm -rf *.out *.o $(EXEC)
