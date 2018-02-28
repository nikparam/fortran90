SUBROUTINE read_params(switch, file_params, params)

	INTEGER, INTENT(IN) :: switch
	CHARACTER, INTENT(IN) :: file_params*40
	DOUBLE PRECISION, INTENT(OUT) :: params(15)
	INTEGER :: i
	DOUBLE PRECISION :: junk(7)

	OPEN(UNIT = 10, FILE = file_params, STATUS = 'OLD',&
	     FORM = 'FORMATTED', ACTION = 'READ')

	DO i = 1, 15
		READ(10,*) junk(:7)

		SELECT CASE(switch)
			CASE(4)
				params(i) = junk(1)
			CASE(5)
				params(i) = junk(2)
			CASE(6)
				params(i) = junk(3)
			CASE(8)
				params(i) = junk(4)
			CASE(12)
				params(i) = junk(5)
			CASE(13)
				params(i) = junk(6)
			CASE(14)
				params(i) = junk(7)
			CASE DEFAULT
				WRITE(*,*) 'WRONG SWITCH'
		END SELECT
	END DO
	
	CLOSE(UNIT = 10)
	RETURN

END SUBROUTINE read_params


SUBROUTINE potential_energy(x,params,V)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: V
	DOUBLE PRECISION :: A, B, C0, C1, C2, D0, D1, D2, &
			    F1, F2, x1, x2, p1, p2, K, &
			    V1, V2, V3, V4, V5

	A = params(1)
	B = params(2)
	C0 = params(3)
	D0 = params(4)
	C1 = params(5)
	D1 = params(6)
	x1 = params(7)
	C2 = params(8)
	D2 = params(9)
	x2 = params(10)
	F1 = params(11)
	p1 = params(12)
	F2 = params(13)
	p2 = params(14)
	K = params(15)
	
	V1 = A / ( 4.0D0 * DATAN(1.0D0) * x * x + B )
	V2 = C0 * DEXP( -D0 * x * x )
	V3 = C1 * DEXP( -D1 * (x - x1) * (x - x1) )
	V4 = C2*DEXP( -D2 * (x + x2) * (x + x2) )
	V5 = F1 * x**p1 + F2 * x**p2
	
	V = V1 + V2 + V3 + V4 + V5 + K
	RETURN

END SUBROUTINE potential_energy

SUBROUTINE diff_potential_energy(x,params,dV)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: dV
	DOUBLE PRECISION :: A, B, C0, C1, C2, D0, D1, D2, &
			    F1, F2, x1, x2, p1, p2, K, &
			    V1, V2, V3, V4, V51, V52

	A = params(1)
	B = params(2)
	C0 = params(3)
	D0 = params(4)
	C1 = params(5)
	D1 = params(6)
	x1 = params(7)
	C2 = params(8)
	D2 = params(9)
	x2 = params(10)
	F1 = params(11)
	p1 = params(12)
	F2 = params(13)
	p2 = params(14)
	K = params(15)

	V1 = -8.0D0 * DATAN(1.0D0) * A * x / ( ( 4.0D0 * DATAN(1.0D0) * x * x + B) * ( 4.0D0 * DATAN(1.0D0) * x * x + B) )
	V2 = -2.0D0 * C0 * D0 * x * DEXP( -D0 * x * x )
	V3 = -2.0D0 * C1 * D1 * (x - x1) * DEXP( -D1 * (x - x1) * (x - x1) )
	V4 = -2.0D0 * C2 * D2 * (x + x2) * DEXP( -D2 * (x + x2) * (x + x2) )
	IF ( p1 .GE. 1.0D0 ) THEN
		V51 = p1 * F1 * x**(p1 - 1.0D0)
	ELSE
		V51 = 0.0D0
	END IF
	IF ( p2 .GE. 1.0D0 ) THEN
		V52 = p2 * F2 * x**(p2 - 1.0D0)
	ELSE
		V52 = 0.0D0
	END IF

	dV = V1 + V2 + V3 + V4 + V51 + V52
	RETURN

END SUBROUTINE diff_potential_energy

SUBROUTINE diff2_potential_energy(x,params,d2V)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: d2V
	DOUBLE PRECISION :: A, B, C0, C1, C2, D0, D1, D2, &
			    F1, F2, x1, x2, p1, p2, K, &
			    V11, V12, V21, V22, V31, V32, V41, V42, V51, V52

	A = params(1)
	B = params(2)
	C0 = params(3)
	D0 = params(4)
	C1 = params(5)
	D1 = params(6)
	x1 = params(7)
	C2 = params(8)
	D2 = params(9)
	x2 = params(10)
	F1 = params(11)
	p1 = params(12)
	F2 = params(13)
	p2 = params(14)
	K = params(15)

	V11 = -8.0D0 * DATAN(1.0D0) * A / ( ( 4.0D0 * DATAN(1.0D0) * x * x + B) * ( 4.0D0 * DATAN(1.0D0) * x * x + B) )
	V12 = 128.0D0 * DATAN(1.0D0)**2 * A * x**2 / ( ( 4.0D0 * DATAN(1.0D0) * x * x + B)**3 ) 
	V21 = -2.0D0 * C0 * D0 * DEXP( -D0 * x * x )
	V22 = 4.0D0 * C0 * D0**2 * x**2 * DEXP( -D0 * x * x )
	V31 = -2.0D0 * C1 * D1 * DEXP( -D1 * (x - x1) * (x - x1) )
	V32 = 4.0D0 * C1 * D1**2 * (x - x1)**2 * DEXP( -D1 * (x - x1) * (x - x1) )
	V41 = -2.0D0 * C2 * D2 * DEXP( -D2 * (x + x2) * (x + x2) )
	V42 = 4.0D0 * C2 * D2**2 * (x + x2)**2 * DEXP( -D2 * (x + x2) * (x + x2) )
	IF ( p1 .GE. 2.0D0 ) THEN
		V51 = p1 * (p1 - 1.0D0) * F1 * x**(p1 - 2.0D0)
	ELSE
		V51 = 0.0D0
	END IF
	IF ( p2 .GE. 2.0D0 ) THEN
		V52 = p2 * (p2 - 1.0D0) * F2 * x**(p2 - 2.0D0)
	ELSE
		V52 = 0.0D0
	END IF

	d2V = V11 + V12 + V21 + V22 + V31 + V32 + V41 + V42 + V51 + V52
	RETURN

END SUBROUTINE diff2_potential_energy

SUBROUTINE diff3_potential_energy(x,params,d3V)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: d3V
	DOUBLE PRECISION :: A, B, C0, C1, C2, D0, D1, D2, &
			    F1, F2, x1, x2, p1, p2, K, &
			    V111, V121, V122, V211, V221, V222, &
			    V311, V321, V322, V411, V421, V422, &
			    V51, V52

	A = params(1)
	B = params(2)
	C0 = params(3)
	D0 = params(4)
	C1 = params(5)
	D1 = params(6)
	x1 = params(7)
	C2 = params(8)
	D2 = params(9)
	x2 = params(10)
	F1 = params(11)
	p1 = params(12)
	F2 = params(13)
	p2 = params(14)
	K = params(15)

	V111 = 128.0D0 * DATAN(1.0D0) * A * x / ( 4.0D0 * DATAN(1.0D0) * x * x + B)**3
	V121 = 256.0D0 * DATAN(1.0D0)**2 * A * x / ( ( 4.0D0 * DATAN(1.0D0) * x * x + B)**3 ) 
	V122 = -3072.0D0 * DATAN(1.0D0)**3 * A * x**3 / ( ( 4.0D0 * DATAN(1.0D0) * x * x + B)**4 )
	V211 = 4.0D0 * C0 * D0**2 * x * DEXP( -D0 * x * x )
	V221 = 8.0D0 * C0 * D0**2 * x * DEXP( -D0 * x * x )
	V222 = -8.0D0 * C0 * D0**3 * x**3 * DEXP( -D0 * x * x )
	V311 = 4.0D0 * C1 * D1**2 * (x - x1) * DEXP( -D1 * (x - x1) * (x - x1) )
	V321 = 8.0D0 * C1 * D1**2 * (x - x1) * DEXP( -D1 * (x - x1) * (x - x1) )
	V322 = -8.0D0 * C1 * D1**3 * (x - x1)**3 * DEXP( -D1 * (x - x1) * (x - x1))
	V411 = -4.0D0 * C2 * D2**2 * (x + x2) * DEXP( -D2 * (x + x2) * (x + x2) )
	V421 = 8.0D0 * C2 * D2**2 * (x + x2) * DEXP( -D2 * (x + x2) * (x + x2) )
	V422 = -8.0D0 * C2 * D2**3 * (x + x2)**3 * DEXP( -D2 * (x + x2) * (x + x2) )
	IF ( p1 .GE. 3.0D0 ) THEN
		V51 = p1 * (p1 - 1.0D0) * (p1 - 2.0D0) * F1 * x**(p1 - 3.0D0) 
	ELSE
		V51 = 0.0D0
	END IF
	IF ( p2 .GE. 3.0D0 ) THEN
		V52 = p2 * (p2 - 1.0D0) * (p2 - 2.0D0) * F2 * x**(p2 - 3.0D0)
	ELSE
		V52 = 0.0D0
	END IF

	d3V = V111 + V121 + V122 + V211 + V221 + V222 + V311 + &
	      V321 + V322 + V411 + V421 + V422 + V51 + V52
	RETURN 

END SUBROUTINE diff3_potential_energy
