SUBROUTINE coeff_matrix( NEQ, m, params, omega, phase, q, p, npts, x, wts, dE, R)

	INTEGER	     :: i
	INTEGER, INTENT(IN) :: NEQ, npts
	DOUBLE PRECISION, INTENT(IN) :: q(NEQ), p(NEQ), m, params(15), omega(NEQ), phase(NEQ), &
					x(npts), wts(npts), dE
	DOUBLE COMPLEX, INTENT(OUT) :: R(NEQ,NEQ)
	DOUBLE PRECISION :: T, dV(NEQ)
	DOUBLE COMPLEX :: ksi(NEQ), eta(NEQ), dksi(NEQ), deta(NEQ), &
			  S(NEQ,NEQ), SI(NEQ,NEQ), &
			  H(NEQ,NEQ), M1(NEQ,NEQ), z_zdot(NEQ,NEQ), &
			  tau1(NEQ,NEQ), tau2(NEQ,NEQ), L(NEQ,NEQ)

	DO i=1,NEQ
		CALL diff_potential_energy(q(i), params, dV(i))
	END DO

	CALL change_var(NEQ, q, p, ksi, eta, omega, phase)
	CALL overlap(NEQ, ksi, eta, omega, S)
	CALL hamiltonian(NEQ, ksi, eta, omega, params, npts, x, wts, dE, S, H, L, M1)

	dksi(1:NEQ) = (/ ( DCMPLX( omega(i) * p(i) / m, -dV(i) ), i=1,NEQ ) /)
	deta(1:NEQ) = (/ ( DCMPLX( -omega(i) * q(i) * p(i) / m, -p(i) * p(i)/m + q(i) * dV(i) ), i=1,NEQ ) /)

	z_zdot(1:NEQ, 1:NEQ) = M1 * SPREAD(dksi(1:NEQ), 1, NEQ) + S * SPREAD(deta(1:NEQ), 1, NEQ)
	CALL inverse(NEQ, S, SI)

	tau1 = MATMUL( SI, H )
	tau2 = ( 0.0D0, 1.0D0 ) * MATMUL( SI, z_zdot )

	R = -( 0.0D0, 1.0D0 ) *  MATMUL( SI, H - (0.0D0, 1.0D0) * z_zdot )

!	DO i = 1, NEQ
!		WRITE(6,20) tau1(i,:)
!	END DO
!
!	WRITE(*,*)
!
!	DO i = 1, NEQ
!		WRITE(6,20) tau2(i,:)
!	END DO
!
!	DO i = 1, NEQ
!		WRITE(6,20) R(i,:)
!	END DO

20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
