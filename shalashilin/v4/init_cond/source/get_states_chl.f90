SUBROUTINE get_states(NumG, omega, phase, params, q, p, eigen_states, eigen_vectors)

	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p(NumG), omega(NumG), params(15)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	INTEGER :: ok, INFO
	DOUBLE PRECISION :: RWORK(2*NumG)
	DOUBLE COMPLEX :: ksi(NumG), eta(NumG), &
			  S(NumG, NumG), H(NumG, NumG), M1(NumG, NumG), &
			  WORK(2*NumG), H_prime(NumG,NumG), &
			  U_prime(NumG,NumG), DUMMY

	CHARACTER :: fname*40
	LOGICAL :: exists_s

	CALL change_var(NumG, q, p, ksi, eta, omega, phase)
	CALL overlap(NumG, ksi, eta, omega, S)
	CALL hamiltonian(NumG, ksi, eta, omega, params, S, H, M1)

	CALL ZPOTRF( 'L', NumG, S, NumG, INFO )
	DO i = 1, NumG
		DO j = i+1, NumG
			S(i,j) = 0.0
		END DO
	END DO

	INFO = 0
	CALL ZTRTRI( 'L', 'N', NumG, S, NumG, INFO )

	DO i = 1, NumG
		DO j = i+1, NumG
			S(i,j) = 0.0
		END DO
	END DO

	H_prime = MATMUL( S , MATMUL( H, TRANSPOSE( CONJG( S ) ) ) )

	CALL ZGEEV('N', 'V', NumG, H_prime, NumG, eigen_states, DUMMY, 1, U_prime, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	print *, '      hamiltonian'
	IF (ok .EQ. 0) THEN
		DO i = 1, NumG
			WRITE(6,'(F14.6)') DBLE( eigen_states(i) )
		END DO
	ELSE
		WRITE(6,*) 'Something went wrong', ok
	END IF

	eigen_vectors = MATMUL( TRANSPOSE( CONJG( S ) ), U_prime )

END SUBROUTINE get_states
