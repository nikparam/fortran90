SUBROUTINE hamiltonian(NEQ, ksi, eta, omega, params, npts, x, wts, dE, S, H, M1)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NEQ, npts
	INTEGER :: i, j
	DOUBLE PRECISION, INTENT(IN) :: omega(NEQ), params(15), x(npts), wts(npts), dE
	DOUBLE COMPLEX, INTENT(IN) ::  ksi(NEQ), eta(NEQ), S(NEQ,NEQ)
	DOUBLE PRECISION :: a(NEQ,NEQ), left_boundary, right_boundary
	DOUBLE COMPLEX ::  b(NEQ,NEQ), c(NEQ,NEQ), &
			   B1(NEQ,NEQ), B2(NEQ,NEQ), M2(NEQ,NEQ), &
			   X0(NEQ), X1(NEQ), X2(NEQ), V(NEQ,NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: H(NEQ,NEQ), M1(NEQ,NEQ)
	DOUBLE PRECISION :: new_x(npts)

	a(1:NEQ,1:NEQ) = 0.5 * ( SPREAD(omega(1:NEQ),1,NEQ) + SPREAD(omega(1:NEQ),2,NEQ) )
	b(1:NEQ,1:NEQ) = SPREAD(ksi(1:NEQ),1,NEQ) + SPREAD(CONJG(ksi(1:NEQ)),2,NEQ)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	M1 = B1 * S
	M2 = B2 * S

	DO i = 1, NEQ
		DO j = i, NEQ
			new_x = shift( npts, x, a(i,j) )
			V(i,j) = SUM( wts * &
				      CONJG( contr_wave_packet( npts, new_x, ksi(i), eta(i) ) ) * &
				      ( potential_map(npts,new_x,params) + dE ) * &
				      contr_wave_packet( npts, new_x, ksi(j), eta(j) ) ) / DSQRT( a(i,j) )
			V(j,i) = CONJG( V(i,j) )
			new_x = 0.0
		END DO
	END DO

	X0 = 0.5 * ( omega - ksi * ksi )
	X1 = omega * ksi
	X2 = -0.5 * omega * omega

	H(:,:) = SPREAD( X0, 1, NEQ ) * S  + &
	    	 SPREAD( X1, 1, NEQ ) * M1 + &
	    	 SPREAD( X2, 1, NEQ ) * M2 

20 	FORMAT( 100("( ",F14.6,SP,F14.6,"i ) ") )

	H(:,:) = H(:,:) + V(:,:)

	CONTAINS

		FUNCTION contr_wave_packet( npts, x, ksi, eta )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega
			DOUBLE COMPLEX :: contr_wave_packet(npts), ksi, eta

			contr_wave_packet(:) = ZEXP( ksi * x(:) + eta )

		END FUNCTION contr_wave_packet

		FUNCTION potential_map(npts, x, params)

			INTEGER :: npts, i
			DOUBLE PRECISION :: x(npts), params(15)
			DOUBLE PRECISION :: potential_map(npts)

			DO i = 1, npts			
				CALL potential_energy( x(i), params, potential_map(i) )
			END DO

		END FUNCTION potential_map

		FUNCTION shift( N, x, a )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a
			DOUBLE PRECISION :: shift(N)

			shift(:) = x(:) / DSQRT( a )

		END FUNCTION shift


END SUBROUTINE hamiltonian
