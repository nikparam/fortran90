SUBROUTINE FEXD( NEQ, T, Y, YDOT, RPAR, IPAR)

	INTEGER	     :: IPAR

	DOUBLE PRECISION :: T

	DOUBLE COMPLEX :: Y(NEQ), YDOT(NEQ), R(NEQ,NEQ), RPAR(NEQ**2)

	R(1:NEQ,1:NEQ) = RESHAPE( RPAR, (/ NEQ, NEQ /))
	YDOT(1:NEQ) = MATMUL( R, Y )
!	WRITE(*,'( 100(F16.8,SP,F16.8) )') YDOT(:)
	RETURN

END SUBROUTINE FEXD
