SUBROUTINE coeff_matrix( NumG, m, params, omega, lambda, ksi, dV, S, H, M1, R )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: m, params(15), omega(NumG), lambda, dV(NumG)
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), S(NumG,NumG), &
				      H(NumG,NumG), M1(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NumG,NumG)

	INTEGER :: i, j
	DOUBLE COMPLEX :: dksi(NumG), deta(NumG), &
			  z_zdot(NumG,NumG), SI(NumG,NumG), &
			  A(NumG,NumG)

	dksi(1:NumG) = (/ ( DCMPLX( omega(i) * DIMAG(ksi(i)), -dV(i) ), i=1,NumG ) /)
	deta(1:NumG) = (/ ( DCMPLX( -DBLE( ksi(i)/m ) * DIMAG(ksi(i)), 0.0 ), i=1,NumG ) /)

!	CALL block_inverse(NumG, lambda, S, SI)
	CALL inverse(NumG, lambda, S, SI)


	z_zdot(:,:) = M1 * SPREAD( dksi, 1, NumG ) + S * SPREAD( deta, 1, NumG )

	A = -( 0.0D0, 1.0D0 ) * ( H - (0.0D0, 1.0D0) * z_zdot )
	R = MATMUL( SI, A )

20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
