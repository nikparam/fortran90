SUBROUTINE potential_energy(x, params, V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: V
	DOUBLE PRECISION :: D_e, r_e, a
	D_e = 1.0D0
	r_e = 0.0D0
	a = 1.0D0 / DSQRT(2.0D0)

	V = D_e * ( 1 - DEXP( -a * ( x - r_e ) ) )**2
	RETURN

END SUBROUTINE potential_energy

SUBROUTINE diff_potential_energy(x, params, dV)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: dV
	DOUBLE PRECISION :: D_e, r_e, a
	D_e = 1.0D0
	r_e = 0.0D0
	a = 1.0D0 / DSQRT(2.0D0)

	dV = 2.0D0 * a * D_e * ( 1 - DEXP( -a * ( x - r_e ) ) ) * &
				     DEXP( -a * ( x - r_e ) )
	RETURN

END SUBROUTINE diff_potential_energy

SUBROUTINE diff2_potential_energy(x, params, d2V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: d2V
	DOUBLE PRECISION :: D_e, r_e, a
	D_e = 1.0D0
	r_e = 0.0D0
	a = 1.0D0 / DSQRT(2.0D0)

	d2V = 2.0D0 * a**2 * D_e * ( 2.0D0 * DEXP( -2.0D0 * a * ( x - r_e ) ) - &
		   	     	  	     DEXP( -a * ( x - r_e ) ))
	RETURN

END SUBROUTINE diff2_potential_energy

SUBROUTINE diff3_potential_energy(x, params, d2V)

	DOUBLE PRECISION, INTENT(IN) :: params(15), x
	DOUBLE PRECISION, INTENT(OUT) :: d2V
	DOUBLE PRECISION :: D_e, r_e, a
	D_e = 1.0D0
	r_e = 0.0D0
	a = 1.0D0 / DSQRT(2.0D0)

	d2V = 2.0D0 * a**3 * D_e * ( -4.0D0 * DEXP( -2.0D0 * a * ( x - r_e ) ) + &
			   	     	      DEXP( -a * ( x - r_e ) ))
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


