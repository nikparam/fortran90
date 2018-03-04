SUBROUTINE hamiltonian(NEQ, ksi, eta, omega, params, S, H, M1)

	INTEGER, INTENT(IN) :: NEQ
	DOUBLE PRECISION, INTENT(IN) :: omega(NEQ), params(15)
	DOUBLE COMPLEX, INTENT(IN) ::  ksi(NEQ), eta(NEQ), S(NEQ,NEQ)
	DOUBLE PRECISION :: a(NEQ,NEQ), V(NEQ), dV(NEQ), d2V(NEQ), d3V(NEQ) 
	DOUBLE COMPLEX ::  b(NEQ,NEQ), c(NEQ,NEQ), &
			   B1(NEQ,NEQ), B2(NEQ,NEQ), B3(NEQ,NEQ), &
			   M2(NEQ,NEQ), M3(NEQ,NEQ), &
			   X0(NEQ), X1(NEQ), X2(NEQ), X3(NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: H(NEQ,NEQ), M1(NEQ,NEQ)

	a(1:NEQ,1:NEQ) = 0.5 * ( SPREAD(omega(1:NEQ),1,NEQ) + SPREAD(omega(1:NEQ),2,NEQ) )
	b(1:NEQ,1:NEQ) = SPREAD(ksi(1:NEQ),1,NEQ) + SPREAD(CONJG(ksi(1:NEQ)),2,NEQ)
	c(1:NEQ,1:NEQ) = SPREAD(eta(1:NEQ),1,NEQ) + SPREAD(CONJG(eta(1:NEQ)),2,NEQ)

!	WRITE(*,*) ksi
!	CALL overlap(NEQ, ksi, eta, a, S)

	B1 = 0.5D0 * b / a
	B2 = 0.25D0 * b**2 / a**2 + 0.5D0 / a
	B3 = 0.125D0 * b**3 / a**3 + 0.75D0 / a**2
	M1 = B1 * S
	M2 = B2 * S
	M3 = B3 * S

!	WRITE(*,*) M1, M2

	DO i = 1, NEQ
		CALL potential_energy( DBLE(ksi(i)) / omega(i), params, V(i) )
		CALL diff_potential_energy( DBLE(ksi(i)) / omega(i), params, dV(i) )
		CALL diff2_potential_energy( DBLE(ksi(i)) / omega(i), params, d2V(i) )
		CALL diff3_potential_energy( DBLE(ksi(i)) / omega(i), params, d3V(i) )
	END DO

!	WRITE(*,*) V, dV, d2V

	X0(1:NEQ) = 0.5D0 * ( omega(1:NEQ) - ksi(1:NEQ) * ksi(1:NEQ) ) + V(1:NEQ) - &
		    dV(1:NEQ) * DBLE(ksi(1:NEQ)) / omega(1:NEQ)  + &
		    0.5D0 * d2V(1:NEQ) * DBLE(ksi(1:NEQ))**2 / omega(1:NEQ)**2 - &
		    d3V(1:NEQ) * DBLE(ksi(1:NEQ))**3 / ( 6.0D0 * omega(1:NEQ)**3 )
	X1(1:NEQ) = omega(1:NEQ) * ksi(1:NEQ) + dV(1:NEQ) - &
		    d2V(1:NEQ) * DBLE(ksi(1:NEQ)) / omega(1:NEQ) + &
		    0.5D0 * d3V(1:NEQ) * DBLE(ksi(1:NEQ))**2 / omega(1:NEQ)**2
	X2(1:NEQ) = -0.5D0 * omega(1:NEQ)**2 + 0.5D0 * d2V(1:NEQ) - &
		    0.5D0 * d3V(1:NEQ) * DBLE(ksi(1:NEQ)) /  omega(1:NEQ)
	X3(1:NEQ) = d3V / 6.0D0

!	DO i = 1, NEQ
!		WRITE(*,*) V(i), dV(i), d2V(i), d3V(i)
!	END DO

	H(1:NEQ,1:NEQ) = SPREAD(X0(1:NEQ),1,NEQ) * S + SPREAD(X1(1:NEQ),1,NEQ) * M1 + &
			 SPREAD(X2(1:NEQ),1,NEQ) * M2 + SPREAD(X3(1:NEQ),1,NEQ) * M3

!	WRITE(*,*) H
!	WRITE(*,*)
	RETURN

END SUBROUTINE hamiltonian
