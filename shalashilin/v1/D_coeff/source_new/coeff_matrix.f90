SUBROUTINE coeff_matrix( NumG, m, params, omega, ksi, S, H, M1, R )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: m, params(15), omega(NumG)
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), S(NumG,NumG), &
				      H(NumG,NumG), M1(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NumG,NumG)

	INTEGER :: i, j
	DOUBLE PRECISION :: T, dV(NumG)
	DOUBLE COMPLEX :: dksi(NumG), deta(NumG), &
			  SI(NumG,NumG), &
			  z_zdot(NumG,NumG)

	DO i=1,NumG
		CALL diff_potential_energy( DBLE( ksi(i) ) / omega(i) , params, dV(i) )
	END DO

	dksi(1:NumG) = (/ ( DCMPLX( omega(i) * DIMAG( ksi(i) ) / m, -dV(i) ), i=1,NumG ) /)
	deta(1:NumG) = (/ ( DCMPLX( -DBLE( ksi(i) ) * DIMAG( ksi(i) ), 0.0 ), i=1,NumG ) /)

	CALL inverse(NumG, S, SI)

	z_zdot(:,:) = M1 * SPREAD( dksi, 1, NumG ) + S * SPREAD( deta, 1, NumG )

	R = -( 0.0D0, 1.0D0 ) *  MATMUL( SI, H - (0.0D0, 1.0D0) * z_zdot )

20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
