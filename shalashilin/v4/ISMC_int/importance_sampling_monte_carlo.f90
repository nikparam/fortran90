SUBROUTINE ISMC_integration( NumG, params, ksi, eta, omega, phase, integral )

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: omega(NumG), phase(NumG), sigma(NumG), params(15)
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), eta(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: integral(NumG,NumG)
	DOUBLE PRECISION :: x1, x2

	CALL change_var(NumG,q,p,ksi,eta,omega,phase)

	sigma = 1.0D0 / DSQRT( 2.0D0 * omega )

	integral = 0.0D0

	DO i = 1, NumG
		DO 10 j=1, 1000
			CALL rand_norm(q(i), p(i), sigma, x1, x2)

			CALL potential_energy(x1, params, V)
			integral(i) = integral + V * DEXP( -omega(i) * (x1 - q(i)) )

			CALL potential_energy(x2, params, V)
			integral(i) = integral + V * DEXP( -omega(i) * (x2 - q(i)) )
10		j = j + 2
	END DO

	integral = integral * 4.0D-3 * sigma 
	
	RETURN

END SUBROUTINE
