SUBROUTINE get_states(NumG, omega, phase, params, q, p, eigen_states, eigen_vectors)

	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p(NumG), omega(NumG), params(15)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_states(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: eigen_vectors(NumG, NumG)

	INTEGER :: ok
	DOUBLE COMPLEX :: WORK( 2 * NumG ), DUMMY
	DOUBLE PRECISION :: RWORK( 2 * NumG )
	DOUBLE COMPLEX :: ksi(NumG), eta(NumG), &
			  S(NumG, NumG), SI(NumG,NumG), &
			  H(NumG, NumG), M1(NumG, NumG), &
			  H_prime(NumG,NumG)

	CHARACTER :: fname*40
	LOGICAL :: exists_s

	CALL change_var(NumG, q, p, ksi, eta, omega, phase)

	CALL overlap(NumG, ksi, eta, omega, S)

	CALL hamiltonian(NumG, ksi, eta, omega, params, S, H, M1)

!	DO i = 1, NumG
!		WRITE(*,20) M1(i,:)
!	END DO
!	WRITE(*,*)
!	DO i = 1, NumG
!		WRITE(*,20) H(i,:)
!	END DO
20 	FORMAT( 100("( ",F14.6,SP,F14.6,"i ) ") )
	CALL inverse(NumG, S, SI)

	H_prime = MATMUL( SI, H )

	CALL ZGEEV('N', 'V', NumG, H_prime, NumG, eigen_states, DUMMY, 1, eigen_vectors, NumG, WORK, 2*NumG, RWORK, 2*NumG, ok)

	print *, '  hamiltonian'
	IF (ok .EQ. 0) THEN
		DO i = 1, NumG
			WRITE(6,'(F14.6)') DBLE( eigen_states(i) )
		END DO
	ELSE
		WRITE(6,*) 'Something went wrong', ok
	END IF

END SUBROUTINE get_states
