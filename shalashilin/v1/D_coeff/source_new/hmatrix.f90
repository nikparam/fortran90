SUBROUTINE hamiltonian(NumG, key, ksi, eta, omega, params, npts, x, wts, dE, S, H, L, M1, M2)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG, key, npts
	DOUBLE PRECISION, INTENT(IN) :: omega(NumG), params(15), dE, x(npts), wts(npts)
	DOUBLE COMPLEX, INTENT(IN) ::  ksi(NumG), eta(NumG), S(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: H(NumG,NumG), L(NumG,NumG), M1(NumG,NumG), M2(NumG,NumG)

	INTEGER :: i, j
	DOUBLE PRECISION :: a(NumG,NumG), V(NumG), dV(NumG), d2V(NumG), d3V(NumG), &
			    left_boundary, right_boundary, new_x(npts)
	DOUBLE COMPLEX ::  b(NumG,NumG), &
			   B1(NumG,NumG), B2(NumG,NumG), B3(NumG,NumG), &
			   M3(NumG,NumG), &
			   X0(NumG), X1(NumG), X2(NumG), X3(NumG), &
			   Vmatrix(NumG,NumG)

	a(1:NumG,1:NumG) = 0.5 * ( SPREAD(omega(1:NumG),1,NumG) + SPREAD(omega(1:NumG),2,NumG) )
	b(1:NumG,1:NumG) = SPREAD(ksi(1:NumG),1,NumG) + SPREAD(CONJG(ksi(1:NumG)),2,NumG)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	M1 = B1 * S
	M2 = B2 * S

	X0 = 0.5 * ( omega - ksi * ksi )
	X1 = omega * ksi
	X2 = -0.5 * omega * omega

	H(:,:) = SPREAD( X0, 1, NumG ) * S  + &
	    	 SPREAD( X1, 1, NumG ) * M1 + &
	    	 SPREAD( X2, 1, NumG ) * M2 

	SELECT CASE( key )
		CASE(0)
			B3 = 0.125D0 * b**3 / a**3 + 0.75D0 * b / a**2
			M3 = B3 * S
			DO i = 1, NumG
				CALL potential_energy( DBLE(ksi(i)) / omega(i), params, V(i) )
				CALL diff_potential_energy( DBLE(ksi(i)) / omega(i), params, dV(i) )
				CALL diff2_potential_energy( DBLE(ksi(i)) / omega(i), params, d2V(i) )
				CALL diff3_potential_energy( DBLE(ksi(i)) / omega(i), params, d3V(i) )
			END DO
		
			X0 = V - dV * DBLE( ksi ) / omega  + &
			     0.5D0 * d2V * DBLE( ksi )**2 / omega**2 - &
			     d3V * DBLE( ksi )**3 / ( 6.0D0 * omega**3 ) + dE
			X1 = dV - d2V * DBLE( ksi ) / omega + &
			     0.5D0 * d3V * DBLE( ksi )**2 / omega**2
			X2 = 0.5D0 * d2V - 0.5D0 * d3V * DBLE( ksi ) / omega
			X3 = d3V / 6.0D0
		
			Vmatrix = SPREAD( X0, 1, NumG ) * S + SPREAD( X1, 1, NumG ) * M1 + &
				  SPREAD( X2, 1, NumG ) * M2 + SPREAD( X3, 1, NumG ) * M3

		CASE(1)
			DO i = 1, NumG
				DO j = i, NumG
					left_boundary = MIN( DBLE(ksi(i))/omega(i) - 6 / DSQRT(omega(i)), &
							     DBLE(ksi(j))/omega(j) - 6 / DSQRT(omega(j)) )
					
					right_boundary = MAX( DBLE(ksi(i))/omega(i) + 6 / DSQRT(omega(i)), &
							      DBLE(ksi(j))/omega(j) + 6 / DSQRT(omega(j)) )
			
					new_x = shift_l( npts, x, left_boundary, right_boundary )

					Vmatrix(i,j) = 0.5 * ( right_boundary - left_boundary ) * &
						       SUM( wts * &
							    CONJG( wave_packet(npts,new_x,ksi(i),eta(i),omega(i)) ) * &
							    ( potential_map(npts,new_x,params) + dE ) * &
							    wave_packet(npts,new_x,ksi(j),eta(j),omega(j)) )
					Vmatrix(j,i) = CONJG( Vmatrix(i,j) )

				END DO
			END DO

!		CASE(2)
!			DO i = 1, NumG
!				DO j = i, NumG
!					new_x = shift_h( npts, x, a(i,j) )
!					Vmatrix(i,j) = SUM( wts * &
!				      			    CONJG( contr_wave_packet( npts, new_x, ksi(i), eta(i) ) ) * &
!						            ( potential_map(npts,new_x,params) + dE ) * &
!				      			    contr_wave_packet( npts, new_x, ksi(j), eta(j) ) ) / DSQRT( a(i,j) )
!					Vmatrix(j,i) = CONJG( Vmatrix(i,j) )
!					new_x = 0.0
!				END DO
!			END DO
	END SELECT

	L = H - Vmatrix
	H = H + Vmatrix

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

		FUNCTION shift_l( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift_l(N)

			shift_l(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift_l

!		FUNCTION contr_wave_packet( npts, x, ksi, eta )
!
!			INTEGER :: npts
!			DOUBLE PRECISION :: x(npts), omega
!			DOUBLE COMPLEX :: contr_wave_packet(npts), ksi, eta
!
!			contr_wave_packet(:) = ZEXP( ksi * x(:) + eta )
!
!		END FUNCTION contr_wave_packet
!
!		FUNCTION shift_h( N, x, a )
!
!			INTEGER :: N
!			DOUBLE PRECISION :: x(N), a
!			DOUBLE PRECISION :: shift_h(N)
!
!			shift_h(:) = x(:) / DSQRT( a )
!
!		END FUNCTION shift_h

END SUBROUTINE hamiltonian
