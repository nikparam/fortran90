PROGRAM test

	IMPLICIT NONE
	INTEGER :: i, N
	DOUBLE PRECISION :: integral, a, b
	DOUBLE PRECISION, ALLOCATABLE :: x(:), wts(:)

	WRITE(*,'(a)') 'INPUT LEVEL OF APPROXIMATION: '
	READ(*,*) N

!	WRITE(*,'(a)') 'INPUT BOUNDARIES (LEFT AND RIGHT): '
!	READ(*,*) a, b

	a = -6
	b = 6

	ALLOCATE( x(N), wts(N) )

	CALL p_quadrature_rule( N, x, wts )

	integral = 0.5 * ( b - a ) * SUM( wts * func( N, shift( N, x, a, b ) ) )
	WRITE(*,*) integral

	CONTAINS
		FUNCTION func( N, x )

			INTEGER :: N
			DOUBLE PRECISION :: x(N)
			DOUBLE PRECISION :: func(N)

			func(:) = DEXP( -0.5 * ( x - 1 ) ** 2 ) * &
				  DEXP( -0.5 * ( x + 1 ) ** 2 ) * &
				  x ** 2

		END FUNCTION func

		FUNCTION shift( N, x, a, b )

			INTEGER :: N
			DOUBLE PRECISION :: x(N), a, b
			DOUBLE PRECISION :: shift(N)

			shift(:) = 0.5 * ( b - a ) * x(:) + 0.5 * ( a + b )

		END FUNCTION shift

END PROGRAM test
