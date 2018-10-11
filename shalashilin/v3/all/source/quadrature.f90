SUBROUTINE quadrature( NumG, my_subroutine, ksi, eta, x, wts, npts, m, omega, params, matrix )

	IMPLICIT NONE

	EXTERNAL :: my_subroutine

	INTEGER, INTENT(IN) :: NumG, npts
	DOUBLE PRECISION, INTENT(IN) :: x(npts), wts(npts), m, omega(NumG), params(15)
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), eta(NumG)

	DOUBLE COMPLEX, INTENT(OUT) :: matrix(NumG,NumG)

	INTEGER :: i, j
	DOUBLE PRECISION :: left_boundary, right_boundary, new_x(npts)


	DO i = 1, NumG
		DO j = i, NumG
			left_boundary = MIN( (DBLE(ksi(i))/omega(i)/m) - 6 / DSQRT(m*omega(i)), &
					     (DBLE(ksi(j))/omega(j)/m) - 6 / DSQRT(m*omega(j)) )
			
			right_boundary = MAX( (DBLE(ksi(i))/omega(i)/m) + 6 / DSQRT(m*omega(i)), &
					      (DBLE(ksi(j))/omega(j)/m) + 6 / DSQRT(m*omega(j)) )
	
			new_x = shift_l( npts, x, left_boundary, right_boundary )

			matrix(i,j) = 0.5 * ( right_boundary - left_boundary ) * &
					      SUM( wts * &
						   CONJG( wave_packet(npts,new_x,ksi(i),eta(i),m,omega(i)) ) * &
						   map_sub(my_subroutine,npts,new_x,params) * &
						   wave_packet(npts,new_x,ksi(j),eta(j),m,omega(j)) )
			IF ( i .NE. j ) THEN

				matrix(j,i) = CONJG( matrix(i,j) )

			END IF

		END DO
	END DO

	RETURN

	CONTAINS
		FUNCTION wave_packet( npts, x, ksi, eta, m, omega )

			INTEGER :: npts
			DOUBLE PRECISION :: x(npts), omega, m
			DOUBLE COMPLEX :: wave_packet(npts), ksi, eta

			wave_packet(:) = ZEXP( -0.5 * m * omega * x(:) * x(:) + ksi * x(:) + eta )

		END FUNCTION wave_packet

		FUNCTION map_sub(sub, npts, x, params)

			EXTERNAL :: sub
			INTEGER :: npts, i
			DOUBLE PRECISION :: x(npts), params(15)
			DOUBLE PRECISION :: map_sub(npts)

			DO i = 1, npts			
				CALL sub( x(i), params, map_sub(i) )
			END DO

		END FUNCTION map_sub

		FUNCTION shift_l( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift_l(N)

			shift_l(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift_l

END SUBROUTINE quadrature
