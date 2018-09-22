SUBROUTINE normalize( NumG, ksi, eta, omega, D )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: omega(NumG)
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), eta(NumG)
	DOUBLE COMPLEX, INTENT(IN OUT) :: D(NumG)
	INTEGER :: i
	DOUBLE PRECISION :: junk(3), N
	DOUBLE COMPLEX :: S(NumG,NumG)

	CALL overlap(NumG,ksi,eta,omega,S)
	
	N = DOT_PRODUCT(D, MATMUL(S,D))
	WRITE(*,*) N
	D(1:NumG) = D(1:NumG) / DSQRT(N)
	RETURN

END SUBROUTINE normalize
