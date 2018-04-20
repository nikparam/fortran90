MODULE GRID

CONTAINS
	SUBROUTINE get_grid( omega, m, factor, params, q_lim, V_lim, NumV, q_array, p_array, omega_array, NumG )

		USE APPEND_ARRAY

		DOUBLE PRECISION, INTENT(IN) :: omega, m, factor, q_lim
		DOUBLE PRECISION, INTENT(IN) :: V_lim(NumV), params(15)
		DOUBLE PRECISION, ALLOCATABLE, INTENT(OUT) :: q_array(:), p_array(:), omega_array(:)
		INTEGER, INTENT(OUT) :: NumG
		INTEGER :: astatus, i, j
		DOUBLE PRECISION :: a, b, p_0, q_0, q, p, V, E, func, corr, divisor_a, divisor_b

!		a = DSQRT(8.0D0 * DATAN(1.0D0) * factor / omega )
!		b = 2.0D0 * DSQRT( 2.0D0 * DATAN(1.0D0) * omega * factor )
!		WRITE(*,*) a, b
		a = 0.20 *  omega
		q = ABS( q_lim )
		q_0 = 0.0D0

		ALLOCATE( q_array(0), p_array(0), omega_array(0), STAT = astatus)

		NumG = 0

15		IF ( q_0 .LT. q ) THEN
			divisor_a = 1.0

			p = DSQRT( 2.0D0 * m * V_lim(NumV) )
!			b = 2.0D0 * DSQRT( 2.0D0 * DATAN(1.0D0) * omega * factor )
			b = 0.20 * omega 
			p_0 = 0.0D0 

!			WRITE(*,*) a, " : "

25			IF ( p_0 .LT. p ) THEN
				divisor_b = 1.0
			
				CALL potential_energy(q_0, params, V)
				E = 0.5D0 * p_0**2 / m + V 

				func = 1.0D0
				DO k = 1, NumV
					func = func * ( E - V_lim(k) )
				END DO
!				WRITE(*,*) "-->" , b, " => ", q_0, p_0, E

				IF ( func .LT. 0.0D0 .AND. ( p_0 .EQ. 0.0D0 .OR. q_0 .EQ. 0.0D0 ) ) THEN
					divisor_b = 1.0
					IF ( E /= 0.0 ) THEN
						divisor_a = 1.0
						IF (q_0 .EQ. 0.0D0) THEN
							DO i = 1, 2
								CALL append(q_array, q_0)
								CALL append(p_array, (-1)**i*p_0)
								CALL append(omega_array, omega)
							END DO
							NumG = NumG + 2
						ELSE IF (p_0 .EQ. 0.0D0) THEN
							DO i = 1, 2
								CALL append(q_array, (-1)**i*q_0)
								CALL append(p_array, p_0)
								CALL append(omega_array, omega)
							END DO
							NumG = NumG + 2
						ELSE
							DO i = 1, 2
								DO j = 1, 2
									CALL append(q_array, (-1)**i*q_0)
									CALL append(p_array, (-1)**j*p_0)
									CALL append(omega_array, omega)
								END DO
							END DO
							NumG = NumG + 4
						END IF
					END IF
				END IF
				p_0 = p_0 + b
				b = b * divisor_b
				GOTO 25
			END IF
			q_0 = q_0 + a
			a = a * divisor_a
			GOTO 15
		END IF

	END SUBROUTINE get_grid

END MODULE GRID


