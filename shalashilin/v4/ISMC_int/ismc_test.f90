PROGRAM ismc_test

	INTEGER :: i, N
	DOUBLE PRECISION :: x1, x2, sigma,eps, params(15), V(2), q(10), V1(10), &
			    integral_in, integral_out, &
			    start_time, end_time
	CHARACTER :: file_params*40

	file_params = '../potential/params.txt'
	CALL read_params(13, file_params, params)

	EPS = 1.0D-9
	q = (/ (-300 + 60 * (i - 1), i = 1, 10)  /)
	sigma = 1 / DSQRT(1.374D-3)
	CALL CPU_TIME(start_time)
	DO i = 1, 10
		integal_in = 0.0D0
		CALL RANDOM_NUMBER(integral_out)
		N = 0
		DO
			CALL rand_norm(q(i), sigma, x1, x2)
			CALL potential_energy(x1, params, V(1))
			CALL potential_energy(x2, params, V(2))

			integral_in = integral_in + V(1) + V(2)

			IF ( ABS(integral_in - integral_out) / (N + 2)  .GT. eps ) THEN
				N = N + 2
				integral_out = integral_in
			ELSE
				EXIT
			END IF

		END DO
		integral_in = integral_in / N
		WRITE(6,10) q(i), integral_in, N
10		FORMAT(F12.6,' ',F12.8,' ',I10.2)
	END DO
	CALL CPU_TIME(end_time)


	WRITE(*,*) (end_time - start_time)

END PROGRAM ismc_test
