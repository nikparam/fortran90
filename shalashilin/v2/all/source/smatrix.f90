SUBROUTINE overlap(NumG, ksi, eta, omega, S)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: omega(NumG)
	DOUBLE COMPLEX, INTENT (IN)  :: ksi(NumG), eta(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: S(NumG,NumG)
	INTEGER :: i, j
	DOUBLE PRECISION :: a(NumG,NumG)
	DOUBLE COMPLEX :: b(NumG,NumG), c(NumG,NumG), e(NumG,NumG)

	a(1:NumG,1:NumG) = 0.5D0 * ( SPREAD(omega(1:NumG),1,NumG) + SPREAD(omega(1:NumG),2,NumG) )
	b(1:NumG,1:NumG) = SPREAD(ksi(1:NumG),1,NumG) + SPREAD(CONJG(ksi(1:NumG)),2,NumG)
	c(1:NumG,1:NumG) = SPREAD(eta(1:NumG),1,NumG) + SPREAD(CONJG(eta(1:NumG)),2,NumG)
	e(1:NumG,1:NumG) = 0.25D0 * b**2 / a + c

	S(1:NumG,1:NumG) = ZEXP( e ) * DSQRT( 4.0D0 * DATAN(1.0D0) / a )

!	DO i = 1, NumG
!	        DO j = 1, NumG
!	       		S(j,i) = ZEXP( e(j,i) ) * DSQRT( 4.0D0 * DATAN(1.0D0) / a(j,i) )
!	        END DO
!	END DO

	RETURN
END SUBROUTINE overlap
