PROGRAM rand_int
	IMPLICIT NONE
	INTEGER, PARAMETER :: N = 10
	INTEGER, PARAMETER :: size_count = 1000
	INTEGER :: i, step_count
	INTEGER,DIMENSION(2) :: x,y,a,b
	INTEGER, DIMENSION(1:size_count) :: count_array
	CHARACTER(LEN = 1) :: string
	
	READ(*,*) string
	DO i = 1, size_count
		x = (/ randint(N), randint(N) /)
		y = (/ randint(N), randint(N) /)
!		WRITE(*,*) 'INITIAL CONDITIONS: (', x(1),',', y(1), '); (', x(2), ',', y(2), ')'
		step_count = 0
		DO
			IF ( string == 'c' ) THEN
				CALL cwandering(N,x(1),y(1),a(1),b(1))
				CALL cwandering(N,x(2),y(2),a(2),b(2))
			ELSE IF ( string == 'w' ) THEN
				CALL wandering(N,a(1),b(1))
				CALL wandering(N,a(2),b(2))
			ELSE
				WRITE(*,*) 'ILLIGAL INPUT'
			END IF
			x = (/ a(1), a(2) /)
			y = (/ b(1), b(2) /)
!			WRITE(*,*) x, y
			IF ( (x(1) ==  x(2)) .AND. (y(1) == y(2)) ) THEN
!				WRITE(*,*) step_count
				EXIT
			END IF
			step_count = step_count + 1
		END DO
		count_array(i) = step_count
	END DO

CONTAINS
	INTEGER FUNCTION randint(N)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: N
		REAL :: r
		CALL RANDOM_NUMBER(r)
		randint = FLOOR( ( N + 1 ) * r )
	END FUNCTION randint

	SUBROUTINE wandering(N,a,b)
		IMPLICIT NONE
		INTEGER, INTENT( IN ) :: N
		INTEGER, INTENT( OUT ) :: a, b

		a = randint(N)
		b = randint(N)
	END SUBROUTINE wandering

	SUBROUTINE cwandering(N,x,y,a,b)
		IMPLICIT NONE
		INTEGER :: step
		INTEGER, INTENT( IN ) :: N, x, y
		INTEGER, INTENT( OUT ) :: a, b
		step = randint(12)

		SELECT CASE ( MOD(step,4) )

			CASE ( 0 )
				IF ( x > 0 ) THEN
					a = x - 1
					b = y
				END IF

			CASE ( 1 )
				IF ( y > 0 ) THEN
					a = x
					b = y - 1
				END IF

			CASE ( 2 )
				IF ( x < 10 ) THEN
					a = x + 1
					b = y
				END IF

			CASE ( 3 )
				IF ( y < 10 ) THEN
					a = x
					b = y + 1
				END IF

			CASE DEFAULT
				WRITE(*,*) 'ERROR in step'

		END SELECT
	END SUBROUTINE cwandering
END PROGRAM rand_int
