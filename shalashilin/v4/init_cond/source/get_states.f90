SUBROUTINE get_states(NumG, omega, phase, params, q, p, eigen_states, eigen_vectors)

	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p(NumG), omega(NumG), params(15)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	INTEGER :: ok
	DOUBLE PRECISION :: RWORK(2*NumG)
	DOUBLE COMPLEX :: ksi(NumG), eta(NumG), &
			  S(NumG, NumG), S_mhalf(NumG,NumG), H(NumG, NumG), M1(NumG, NumG), &
			  s_eigen(NumG), TL(NumG, NumG), TR(NumG, NumG), &
			  WORK(2*NumG), H_prime(NumG,NumG), H_dprime(NumG,NumG), &
			  U_prime(NumG,NumG), DUMMY, TEMP(1:NumG,1:NumG)


	CALL change_var(NumG, q, p, ksi, eta, omega, phase)
	CALL overlap(NumG, ksi, eta, omega, S)
	TEMP(1:NumG, 1:NumG) = S(1:NumG, 1:NumG)
	CALL hamiltonian(NumG, ksi, eta, omega, params, S, H, M1)

	ok = 0
	CALL ZGEEV('N', 'V', NumG, S, NumG, s_eigen, DUMMY, 1, TR, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

20	FORMAT(F16.8,' ', F16.8)

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

!	TEMP = MATMUL( TRANSPOSE( CONJG(eigen_vectors) ), MATMUL(TEMP, eigen_vectors) )
!	DO i = 1, NumG
!		WRITE(*,*) TEMP(i,:)
!	END DO

END SUBROUTINE get_states
