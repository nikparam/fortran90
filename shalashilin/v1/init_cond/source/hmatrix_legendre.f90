SUBROUTINE hamiltonian(NEQ, ksi, eta, omega, params, S, H, M1)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NEQ
	INTEGER :: i, j, npts
	DOUBLE PRECISION, INTENT(IN) :: omega(NEQ), params(15)
	DOUBLE COMPLEX, INTENT(IN) ::  ksi(NEQ), eta(NEQ), S(NEQ,NEQ)
	DOUBLE PRECISION :: a(NEQ,NEQ), left_boundary, right_boundary, t_start, t_finish
	DOUBLE COMPLEX ::  b(NEQ,NEQ), c(NEQ,NEQ), &
			   B1(NEQ,NEQ), B2(NEQ,NEQ), M2(NEQ,NEQ), &
			   X0(NEQ), X1(NEQ), X2(NEQ), V(NEQ,NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: H(NEQ,NEQ), M1(NEQ,NEQ)
	DOUBLE PRECISION, ALLOCATABLE :: x(:), new_x(:), wts(:)

	a(1:NEQ,1:NEQ) = 0.5 * ( SPREAD(omega(1:NEQ),1,NEQ) + SPREAD(omega(1:NEQ),2,NEQ) )
	b(1:NEQ,1:NEQ) = SPREAD(ksi(1:NEQ),1,NEQ) + SPREAD(CONJG(ksi(1:NEQ)),2,NEQ)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	M1 = B1 * S
	M2 = B2 * S

	npts = 50
	ALLOCATE( x(npts), new_x(npts), wts(npts) )
	CALL p_quadrature_rule( npts, x, wts )

	DO i = 1, NEQ
		DO j = i, NEQ
			left_boundary = MIN( DBLE(ksi(i))/omega(i) - 6 / DSQRT(omega(i)), &
					     DBLE(ksi(j))/omega(j) - 6 / DSQRT(omega(j))   )
					
			right_boundary = MAX( DBLE(ksi(i))/omega(i) + 6 / DSQRT(omega(i)), &
					      DBLE(ksi(j))/omega(j) + 6 / DSQRT(omega(j))   )

			
			new_x = shift( npts, x, left_boundary, right_boundary )

			V(i,j) = 0.5 * ( right_boundary - left_boundary ) * &
				 SUM( wts * &
				      CONJG( wave_packet(npts,new_x,ksi(i),eta(i),omega(i)) ) * &
				      potential_map(npts,new_x,params) * &
				      wave_packet(npts,new_x,ksi(j),eta(j),omega(j)) )
			V(j,i) = CONJG( V(i,j) )
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

		FUNCTION wave_packet( npts, x, ksi, eta, omega )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega
			DOUBLE COMPLEX :: wave_packet(npts), ksi, eta

			wave_packet(:) = ZEXP( -0.5 * omega * x(:) * x(:) + ksi * x(:) + eta )

		END FUNCTION wave_packet

		FUNCTION potential_map(npts, x, params)

			INTEGER :: npts, i
			DOUBLE PRECISION :: x(npts), params(15)
			DOUBLE PRECISION :: potential_map(npts)

			DO i = 1, npts			
				CALL potential_energy( x(i), params, potential_map(i) )
			END DO

		END FUNCTION potential_map

		FUNCTION shift( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift(N)

			shift(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift


END SUBROUTINE hamiltonian
