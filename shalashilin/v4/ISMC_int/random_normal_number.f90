SUBROUTINE rand_norm( x, sigma, x1, x2 )

	DOUBLE PRECISION, INTENT(IN) :: x, sigma
	DOUBLE PRECISION, INTENT(OUT) :: x1, x2
	DOUBLE PRECISION :: u, v, y1, y2


	CALL RANDOM_NUMBER(u)
	CALL RANDOM_NUMBER(v)

	y1 = DSQRT( -2.0D0 * DLOG(u) ) * DCOS( 8.0D0 * DATAN(1.0D0) * v )
	y2 = DSQRT( -2.0D0 * DLOG(u) ) * DSIN( 8.0D0 * DATAN(1.0D0) * v )

	x1 = x + y1 * sigma
	x2 = x + y2 * sigma

	RETURN

END SUBROUTINE rand_norm
