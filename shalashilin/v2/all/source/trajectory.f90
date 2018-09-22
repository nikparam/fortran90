SUBROUTINE trajectory(NumG, q, p_0, Step, m, q_next )

	IMPLICIT NONE

	EXTERNAL :: FEXq

	INTEGER, INTENT(IN) :: NumG 
	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p_0, Step, m
	DOUBLE PRECISION, INTENT(OUT) :: q_next(NumG)

	INTEGER :: itraj, ITASK, ISTATE, IOPT, LZW, LRW, IWORK(50), &
		   LIW, MF, IPAR, ITOL 
	DOUBLE PRECISION :: RWORK(100), &
			    DUMMY(1,1), RPAR(2), &
			    ATOL, RTOL 
	DOUBLE PRECISION :: T, TOUT
	DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-5, EPS_step = 1.0D-4
	DOUBLE COMPLEX :: ZWORK(100), Y

	ITOL = 1
	RTOL = 1.0D-14
	ATOL = 1.0D-14
	IOPT = 0
	LZW = 100
	LRW = 100
	LIW = 50
	MF = 20

	RPAR(1) = p_0
	RPAR(2) = m
	IPAR = NumG

	DO itraj = 1, NumG

		T = 0.0D0
		TOUT = Step

		ISTATE = 1
		ITASK = 1
		ZWORK = ( 0.0D0, 0.0D0 )
		RWORK = 0.0D0
		IWORK = 0

		q_next(itraj) = q(itraj)
		Y = DCMPLX( q_next(itraj), 0.0D0 )

		CALL ZVODE(FEXq, 1, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, &
			   IOPT, ZWORK, LZW, RWORK, LRW, IWORK, LIW, DUMMY, MF, RPAR, IPAR)

		q_next(itraj) = DBLE( Y )
		Y = 0.0

		IF ( ISTATE .LT. 0 ) GO TO 15
	END DO		
	RETURN

	STOP
15	WRITE(6,16) ISTATE
16	FORMAT(///'  Error halt: ISTATE =',I3)
	STOP


END SUBROUTINE trajectory
