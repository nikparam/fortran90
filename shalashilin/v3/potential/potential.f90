PROGRAM potential

	IMPLICIT NONE
	INTEGER :: switch
	DOUBLE PRECISION, PARAMETER :: step = 0.1D0
	DOUBLE PRECISION :: x, V, V1, &
			    dV, dV1, dVn, d2V, dV2, d2Vn, d3V, d3Vn
	INTEGER :: i
	DOUBLE PRECISION :: params(15), junk(6)
	CHARACTER :: fname*40

	READ(*,*) switch
	fname = './params.txt'
	CALL read_params(switch, fname, params(:15))

	WRITE(*,*) params
	x = 0.0D0
	DO 40 i = 1, 100000
		CALL potential_energy(x, params, V)
		IF ( i .GT. 1 ) THEN
			dVn = ( V - V1 ) / step
		END IF

		CALL diff_potential_energy(x, params, dV)
		IF ( i .GT. 1 ) THEN
			d2Vn = ( dV - dV1 ) / step
		END IF
		
		CALL diff2_potential_energy(x, params, d2V)
		IF ( i .GT. 1 ) THEN
			d3Vn = ( d2V - dV2 ) / step
		END IF

		CALL diff3_potential_energy(x, params, d3V)
		IF ( i .GT. 1) THEN
			WRITE(6,70) x, V, dV, dVn, d2V, d2Vn, d3V, d3Vn
		ELSE
			WRITE(6,70) x, V, dV, 0.0D0, d2V, 0.0D0, d3V, 0.0D0
		END IF

70		FORMAT( F12.6, 7F12.8 )

		V1 = V
		dV1 = dV
		dV2 = d2V

	40 x = x + step

!	DO i = 1,15
!		WRITE(*,*) params(i)
!	END DO
END PROGRAM potential
