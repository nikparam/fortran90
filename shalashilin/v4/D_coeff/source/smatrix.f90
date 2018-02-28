SUBROUTINE overlap(NEQ, x, y, const, M)

	INTEGER :: i, j
	INTEGER, INTENT(IN) :: NEQ
	DOUBLE PRECISION, INTENT(IN) :: const(NEQ)
	DOUBLE COMPLEX, INTENT (IN)  :: x(NEQ), y(NEQ)
	DOUBLE PRECISION :: a(NEQ,NEQ)
	DOUBLE COMPLEX :: b(NEQ,NEQ), c(NEQ,NEQ), e(NEQ,NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: M(NEQ,NEQ)

	a(1:NEQ,1:NEQ) = 0.5D0 * ( SPREAD(const(1:NEQ),1,NEQ) + SPREAD(const(1:NEQ),2,NEQ) )
	b(1:NEQ,1:NEQ) = SPREAD(x(1:NEQ),1,NEQ) + SPREAD(CONJG(x(1:NEQ)),2,NEQ)
	c(1:NEQ,1:NEQ) = SPREAD(y(1:NEQ),1,NEQ) + SPREAD(CONJG(y(1:NEQ)),2,NEQ)
	e(1:NEQ,1:NEQ) = 0.25D0 * b**2 / a + c


!	DO i = 1, NEQ
!		WRITE(*,*) e(i,:)
!	END DO

	M(1:NEQ,1:NEQ) = (0.0D0, 0.0D0)
	DO i = 1, NEQ
	        DO j = 1, NEQ
	       		M(j,i) = EXP( e(j,i) ) * DSQRT( 4.0D0 * DATAN(1.0D0) / a(j,i) )
	        END DO
	END DO

!	DO i = 1, NEQ
!		WRITE(*,*) M(i,:)
!	END DO

!	M = M * DSQRT( 4.0D0 * DATAN(1.0D0) / a ) 
	RETURN
END SUBROUTINE overlap
