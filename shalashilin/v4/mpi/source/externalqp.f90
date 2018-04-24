SUBROUTINE FEXqp( NEQ, T, Y, YDOT, RPAR, IPAR)

	REAL*8 :: Y(NEQ), YDOT(NEQ), dV, T
	REAL*8 :: RPAR(16) 
!	WRITE(*,*) T, Y
	CALL diff_potential_energy( Y(1), RPAR(2:16), dV ) 
	YDOT(2) = -dV
	YDOT(1) = Y(2) / RPAR(1)
	RETURN
END SUBROUTINE
