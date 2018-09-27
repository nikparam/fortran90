SUBROUTINE trajectories( Step, NumG, q, q_0, p_0, m, params, q_next, p_0_next )

	IMPLICIT NONE

	EXTERNAl :: FEXp

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: Step, m, q(NumG), q_0, p_0, params(15)
	DOUBLE PRECISION, INTENT(OUT) :: q_next(NumG), p_0_next

	INTEGER :: i

	INTEGER :: ITASK, ISTATE, IOPT, LZW, LRW, IWORK(50), LIW, MF, IPAR, ITOL

	DOUBLE PRECISION :: T, TOUT, RWORK(100), DUMMY(1,1), RPAR(16), ATOL, RTOL

	DOUBLE COMPLEX :: ZWORK(100), Y

	CALL trajectory(NumG, q, p_0, Step, m, q_next)

	p_0_next = p_0

	T = 0.0D0
	TOUT = Step

	ITOL = 1
	RTOL = 1.0D-15
	ATOL = 1.0D-15
	IOPT = 0
	LZW = 100
	LRW = 100
	LIW = 50
	MF = 20
	ISTATE = 1
	ITASK = 1
	ISTATE = 1
	ZWORK = ( 0.0D0, 0.0D0 )
	RWORK = 0.0D0
	IWORK = 0

	RPAR(1) = q_0 
	RPAR(2:16) = params(:)

	Y = DCMPLX( p_0_next, 0.0D0 )

	CALL ZVODE(FEXp, 1, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, &
		   IOPT, ZWORK, LZW, RWORK, LRW, IWORK, LIW, DUMMY, MF, RPAR, IPAR)

	p_0_next = DBLE( Y )
	Y = 0.0D0

	RETURN

END SUBROUTINE trajectories

