SUBROUTINE normalize(NumG,q,p,omega,phase,D)

	INTEGER, INTENT(IN) :: NumG
	INTEGER :: i
	DOUBLE PRECISION :: junk(3), N
	DOUBLE PRECISION, INTENT(IN) :: phase(NumG), omega(NumG), q(NumG), p(NumG)
	DOUBLE COMPLEX :: ksi(NumG), eta(NumG), S(NumG,NumG)
	DOUBLE COMPLEX, INTENT(IN OUT) :: D(NumG)

	CALL change_var(NumG,q,p,ksi,eta,omega,phase)
	CALL overlap(NumG,ksi,eta,omega,S)
	
	N = DOT_PRODUCT(D, MATMUL(S,D))
	WRITE(*,*) N
	D(1:NumG) = D(1:NumG) / DSQRT(N)
	RETURN

END SUBROUTINE normalize
