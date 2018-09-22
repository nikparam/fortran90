SUBROUTINE hamiltonian(NEQ, ksi, eta, coeff, omega, params, npts, x, wts, dE, S, H, L, M1, Y)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: NEQ, npts
	INTEGER :: i, j
	DOUBLE PRECISION, INTENT(IN) :: omega(NEQ), params(15), x(npts), wts(npts), dE
	DOUBLE COMPLEX, INTENT(IN) ::  ksi(NEQ), eta(NEQ), S(NEQ,NEQ), coeff(NEQ)
	DOUBLE PRECISION :: a(NEQ,NEQ), left_boundary, right_boundary, t_start, t_finish
	DOUBLE COMPLEX ::  b(NEQ,NEQ), c(NEQ,NEQ), &
			   B1(NEQ,NEQ), B2(NEQ,NEQ), B3(NEQ,NEQ), M2(NEQ,NEQ), M3(NEQ,NEQ), &
			   X0(NEQ), X1(NEQ), X2(NEQ), X3(NEQ),&
			   V(NEQ,NEQ), &
			   t1(NEQ,NEQ), t2(NEQ,NEQ), t3(NEQ,NEQ), t4(NEQ,NEQ), &
			   temp1(NEQ), temp2(NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: H(NEQ,NEQ), L(NEQ,NEQ), M1(NEQ,NEQ), Y(NEQ)
	DOUBLE PRECISION :: new_x(npts)

	a(1:NEQ,1:NEQ) = 0.5 * ( SPREAD(omega(1:NEQ),1,NEQ) + SPREAD(omega(1:NEQ),2,NEQ) )
	b(1:NEQ,1:NEQ) = SPREAD(ksi(1:NEQ),1,NEQ) + SPREAD(CONJG(ksi(1:NEQ)),2,NEQ)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	B3 = 0.125D0 * b**3 / a**3 + 0.75D0 * b / a**2
	M1 = B1 * S
	M2 = B2 * S
	M3 = B3 * S

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
				      ( potential_map(npts,new_x,params) + dE ) * &
				      wave_packet(npts,new_x,ksi(j),eta(j),omega(j)) )
			V(j,i) = CONJG( V(i,j) )
		END DO
	END DO

	X0 = 0.5 * ( omega - ksi * ksi )
	X1 = omega * ksi
	X2 = -0.5 * omega * omega
	X3 = -0.5 * ksi

	H(:,:) = SPREAD( X0, 1, NEQ ) * S  + &
	    	 SPREAD( X1, 1, NEQ ) * M1 + &
	    	 SPREAD( X2, 1, NEQ ) * M2 

20 	FORMAT( 100("( ",F14.6,SP,F14.6,"i ) ") )

	L(:,:) = H(:,:) - V(:,:)
	H(:,:) = H(:,:) + V(:,:)

	t1(:,:) = MATMUL( H, M1 )
	t2(:,:) = MATMUL( M1, H )
	t3(:,:) = SPREAD( X0 + 0.5 * omega, 1, NEQ ) * M1 + &
	    	  SPREAD( X1, 1, NEQ ) * M2 + &
	    	  SPREAD( X2, 1, NEQ ) * M3 + &
		  SPREAD( X3, 1, NEQ ) * S
	t4(:,:) = SPREAD( X0, 1, NEQ ) * M1 + &
	    	  SPREAD( X1, 1, NEQ ) * M2 + &
	    	  SPREAD( X2, 1, NEQ ) * M3 

	temp1(:) = MATMUL( CONJG(coeff), t3 - t1 )
	temp2(:) = MATMUL( t4 - t2, coeff )

	Y(:) = 0.0
	DO i = 1, NEQ
		DO j = 1, NEQ
			Y( i ) = Y( i ) + ( t3(j,i) - t1(j,i) ) * coeff( i ) * conjg( coeff( j ) ) + &
					  ( t4(i,j) - t2(i,j) ) * CONJG( coeff( i ) ) * coeff( j )
		END DO
	END DO

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
