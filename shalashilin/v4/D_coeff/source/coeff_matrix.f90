SUBROUTINE coeff_matrix( NEQ, m, params, omega, phase, q, p, R)

	INTEGER	     :: i 
	INTEGER, INTENT(IN) :: NEQ
	DOUBLE PRECISION, INTENT(IN) :: q(NEQ), p(NEQ), m, params(15), omega(NEQ), phase(NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NEQ,NEQ)
	DOUBLE PRECISION :: T, dV(NEQ)
	DOUBLE COMPLEX :: ksi(NEQ), eta(NEQ), dksi(NEQ), deta(NEQ), &
			  S(NEQ,NEQ), SI(NEQ,NEQ), &
			  H(NEQ,NEQ), M1(NEQ,NEQ), z_zdot(NEQ,NEQ), tau(NEQ,NEQ)


	DO i=1,NEQ
		CALL diff_potential_energy(q(i), params, dV(i))
	END DO

	CALL change_var(NEQ, q, p, ksi, eta, omega, phase)
	CALL overlap(NEQ, ksi, eta, omega, S)
	CALL hamiltonian(NEQ, ksi, eta, omega, params, S, H, M1)

	dksi(1:NEQ) = (/ ( DCMPLX( omega(i) * p(i) / m, -dV(i) ), i=1,NEQ ) /)
	deta(1:NEQ) = (/ ( DCMPLX( -omega(i) * q(i) * p(i) / m, 0.0D0 ), i=1,NEQ ) /)

	z_zdot(1:NEQ, 1:NEQ) = M1 * SPREAD(dksi(1:NEQ), 1, NEQ) + S * SPREAD(deta(1:NEQ), 1, NEQ)

!	tmp = S
	CALL inverse(NEQ, S, SI)

!	tmp = MATMUL(tmp,SI)
!	DO i = 1, NEQ
!		print *, tmp(i, :)
!	END DO

	tau = H - (0.0D0, 1.0D0) * z_zdot

	R = -( 0.0D0, 1.0D0 ) *  MATMUL( SI, tau )

	RETURN
END SUBROUTINE coeff_matrix
