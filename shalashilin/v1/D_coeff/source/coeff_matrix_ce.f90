SUBROUTINE coeff_matrix( NEQ, m, params, omega, phase, q, p, C, npts, x, wts, E, dE, R)

	INTEGER	     :: i
	INTEGER, INTENT(IN) :: NEQ, npts
	DOUBLE PRECISION, INTENT(IN) :: q(NEQ), p(NEQ), m, params(15), omega(NEQ), phase(NEQ), &
					x(npts), wts(npts), E, dE
	DOUBLE COMPLEX, INTENT(IN) :: C(NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NEQ,NEQ)
	DOUBLE PRECISION :: T, dV(NEQ)
	DOUBLE COMPLEX :: ksi(NEQ), eta(NEQ), dksi(NEQ), deta(NEQ), &
			  S(NEQ,NEQ), SI(NEQ,NEQ), &
			  H(NEQ,NEQ), M1(NEQ,NEQ), z_zdot(NEQ,NEQ), &
			  L(NEQ,NEQ), d, Y(NEQ), temp(NEQ,NEQ), corr(NEQ,NEQ)

	DO i=1,NEQ
		CALL diff_potential_energy(q(i), params, dV(i))
	END DO

	CALL change_var(NEQ, q, p, ksi, eta, omega, phase)
	CALL overlap(NEQ, ksi, eta, omega, S)
	CALL hamiltonian(NEQ, ksi, eta, C, omega, params, npts, x, wts, dE, S, H, L, M1, Y)

	dksi(1:NEQ) = (/ ( DCMPLX( omega(i) * p(i) / m, -dV(i) ), i=1,NEQ ) /)
	deta(1:NEQ) = (/ ( DCMPLX( -omega(i) * q(i) * p(i) / m, 0.0 ), i=1,NEQ ) /)

	z_zdot(1:NEQ, 1:NEQ) = M1 * SPREAD(dksi(1:NEQ), 1, NEQ) + S * SPREAD(deta(1:NEQ), 1, NEQ)
	temp(:,:) = S(:,:)
	CALL inverse(NEQ, S, SI)

	d = DOT_PRODUCT( C, MATMUL( H, MATMUL( SI, MATMUL( SI, MATMUL( H, C ) ) ) ) ) - E * E

	corr(:,:) = 2 * DBLE( DOT_PRODUCT( Y, dksi ) ) * ( H(:,:) - E * S(:,:) ) / d

	R = -( 0.0D0, 1.0D0 ) *  MATMUL( SI, H - (0.0D0, 1.0D0) * z_zdot - corr )


20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
