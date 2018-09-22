SUBROUTINE potential_energy(x, params, V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: V
	DOUBLE PRECISION :: r_e
	r_e = 0.0D0

	V = DEXP( - ( x - r_e ) ) 
	RETURN

END SUBROUTINE potential_energy

SUBROUTINE diff_potential_energy(x, params, dV)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: dV
	DOUBLE PRECISION :: r_e
	r_e = 0.0D0

	dV = -DEXP( - ( x - r_e ) ) 
	RETURN

END SUBROUTINE diff_potential_energy

SUBROUTINE diff2_potential_energy(x, params, d2V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: d2V
	DOUBLE PRECISION :: r_e
	r_e = 0.0D0

	d2V = DEXP( - ( x - r_e ) ) 
	RETURN

END SUBROUTINE diff2_potential_energy

SUBROUTINE diff3_potential_energy(x, params, d3V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: d3V
	DOUBLE PRECISION :: r_e
	r_e = 0.0D0

	d3V = -DEXP( - ( x - r_e ) ) 
	RETURN

END SUBROUTINE diff3_potential_energy

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


