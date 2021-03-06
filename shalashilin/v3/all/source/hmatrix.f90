SUBROUTINE hamiltonian(NumG, key, ksi, eta, m, omega, params, npts, x, wts, S, H, L, M1, M2)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG, key, npts
	DOUBLE PRECISION, INTENT(IN) :: m, omega(NumG), params(15), x(npts), wts(npts)
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

	a(1:NumG,1:NumG) = 0.5 * m * ( SPREAD(omega(1:NumG),1,NumG) + SPREAD(omega(1:NumG),2,NumG) )
	b(1:NumG,1:NumG) = SPREAD(ksi(1:NumG),1,NumG) + SPREAD(CONJG(ksi(1:NumG)),2,NumG)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	M1 = B1 * S
	M2 = B2 * S

	X0 = 0.5 * ( omega - ksi * ksi / m )
	X1 = omega * ksi
	X2 = -0.5 * m * omega * omega

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
		
			X0 = V - dV * ( DBLE( ksi ) / omega / m )  + &
			     0.5D0 * d2V * ( DBLE( ksi ) / omega / m )**2 - &
			     d3V * ( DBLE( ksi )**3 /  omega / m )**3 / 6.0D0
			X1 = dV - d2V * ( DBLE( ksi ) / omega / m ) + &
			     0.5D0 * d3V * ( DBLE( ksi ) / omega / m )**2
			X2 = 0.5D0 * d2V - 0.5D0 * d3V * ( DBLE( ksi ) / omega / m )
			X3 = d3V / 6.0D0
		
			Vmatrix = SPREAD( X0, 1, NumG ) * S + SPREAD( X1, 1, NumG ) * M1 + &
				  SPREAD( X2, 1, NumG ) * M2 + SPREAD( X3, 1, NumG ) * M3

		CASE(1)

			CALL quadrature(NumG, potential_energy, ksi, eta, x, wts, npts, m, omega, params, Vmatrix )

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
