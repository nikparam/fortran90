! Here are external subroutines:
! vector of derivatives
	SUBROUTINE FEX(NEQ, T, Y, YDOT)
		IMPLICIT NONE
		INTEGER :: NEQ
		DOUBLE PRECISION :: T, Y(3), YDOT(3)
		YDOT(1) = -0.04D0 * Y(1) + 1.D4 * Y(2) * Y(3)
		YDOT(3) = 3.D7 * Y(2) * Y(2)
		YDOT(2) = - YDOT(1) - YDOT(3)
		RETURN
	END SUBROUTINE FEX

! Jacobian matrix
	SUBROUTINE JEX(NEQ, T, Y, ML, MU, PD, NRPD)
		IMPLICIT NONE
		INTEGER :: NEQ, ML, MU, NRPD
		DOUBLE PRECISION :: PD(NRPD, 3), T, Y(3)
		PD(1,1) = -0.04D0
		PD(1,2) = 1.D4*Y(3)
		PD(1,3) = 1.D4*Y(2)
		PD(2,1) = 0.04D0
		PD(2,3) = -PD(1,3)
		PD(3,2) = 6.D7*Y(2)
		PD(2,2) = -PD(1,2) - PD(3,2)
		RETURN
	END SUBROUTINE JEX

