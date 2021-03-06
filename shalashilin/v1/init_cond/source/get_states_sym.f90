SUBROUTINE get_states(NumG, omega, phase, params, q, p, eigen_states, eigen_vectors)

	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p(NumG), omega(NumG), params(15)
	DOUBLE PRECISION :: V_junk
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	INTEGER :: ok
	DOUBLE PRECISION :: RWORK(2*NumG)
	DOUBLE COMPLEX :: ksi(NumG), eta(NumG), &
			  S(NumG, NumG), S_mhalf(NumG,NumG), H(NumG, NumG), M1(NumG, NumG), &
			  s_eigen(NumG), TL(NumG, NumG), TR(NumG, NumG), &
			  WORK(2*NumG), H_prime(NumG,NumG), H_dprime(NumG,NumG), &
			  U_prime(NumG,NumG), DUMMY, TEMP(1:NumG,1:NumG)

	CHARACTER :: fname*40
	LOGICAL :: exists_s

	CALL change_var(NumG, q, p, ksi, eta, omega, phase)
	CALL overlap(NumG, ksi, eta, omega, S)
	TEMP(1:NumG, 1:NumG) = S(1:NumG, 1:NumG)
	CALL hamiltonian(NumG, ksi, eta, omega, params, S, H, M1)

	ok = 0
	CALL ZGEEV('N', 'V', NumG, S, NumG, s_eigen, DUMMY, 1, TR, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

20	FORMAT(100(F16.8,' ', F16.8))

	fname = 'lowdin.out'
	INQUIRE(FILE = fname, EXIST = exists_s)
	IF ( exists_s ) THEN
		OPEN(UNIT = 20, FILE = fname, FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 20, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 20, FILE = fname, FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	WRITE(20,*) "	q		p		E"
	DO i = 1, NumG
		CALL potential_energy( q(i), params, V_junk )
		WRITE(20,'(F14.6," ",F14.6," ",F14.6)') q(i), p(i), 0.5D0 * p(i)**2 / m + V_junk
	END DO

	WRITE(20,*)

	DO i = 1, NumG
		WRITE(20,*) '___________________________________________'
		DO j = 1, NumG
			WRITE(20,'("U[",I2.2,",",I2.2,"]= ",F14.6,SP,F14.6,"i")') j, i, TR(j,i)
		END DO
		WRITE(20,*)
	END DO

	S_mhalf(1:NumG,1:NumG) = (0.0D0, 0.0D0)
	DO i = 1, NumG
		S_mhalf(i,i) = 1.0D0 / DSQRT( DBLE( s_eigen(i) ) )
	END DO

	H_prime = MATMUL( TRANSPOSE( CONJG( TR ) ), MATMUL( H, TR ) )
	H_dprime = MATMUL( S_mhalf, MATMUL( H_prime, S_mhalf ) )

	WORK = ( 0.0D0, 0.0D0 )
	RWORK = 0.0D0
	DUMMY = ( 0.0D0, 0.0D0 )
	ok = 0

	CALL ZGEEV('N', 'V', NumG, H_dprime, NumG, eigen_states, DUMMY, 1, U_prime, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	print *, '      overlap        hamiltonian'
	IF (ok .EQ. 0) THEN
		DO i = 1, NumG
			WRITE(6,20) DBLE( s_eigen(i) ), DBLE( eigen_states(i) )
		END DO
	ELSE
		WRITE(6,*) 'Something went wrong', ok
	END IF

	eigen_vectors = MATMUL( TR, MATMUL( S_mhalf, U_prime ) )

END SUBROUTINE get_states
