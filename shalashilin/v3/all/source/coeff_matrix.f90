SUBROUTINE coeff_matrix( NumG, mean, m, params, omega, lambda, ksi, q_0, S, H, M1, R )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG, mean

	DOUBLE PRECISION, INTENT(IN) :: m, params(15), omega(NumG), q_0, lambda

	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), S(NumG,NumG), &
				      H(NumG,NumG), M1(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NumG,NumG)

	INTEGER :: i, j
	DOUBLE PRECISION :: T, dV(NumG)
	DOUBLE COMPLEX :: dksi(NumG), deta(NumG), &
			  z_zdot(NumG,NumG), SI(NumG,NumG), &
			  A(NumG,NumG)

	IF ( mean .EQ. 0 ) THEN
		DO i = 1, NumG
			CALL diff_potential_energy( (DBLE(ksi(i))/omega(i)/m), params, dV(i) ) 
		END DO

	ELSE
		DO i = 1, NumG
			CALL diff_potential_energy( q_0, params, dV(i) )
		END DO
	END IF

	dksi(1:NumG) = (/ ( DCMPLX( omega(i) * DIMAG(ksi(i)), -dV(i) ), i=1,NumG ) /)
	deta(1:NumG) = (/ ( DCMPLX( -DBLE( ksi(i)/m ) * DIMAG(ksi(i)), 0.0 ), i=1,NumG ) /)
	CALL inverse(NumG, lambda, S, SI)

	z_zdot(:,:) = M1 * SPREAD( dksi, 1, NumG ) + S * SPREAD( deta, 1, NumG )

	A = -( 0.0D0, 1.0D0 ) * ( H - (0.0D0, 1.0D0) * z_zdot )
	R = MATMUL( SI, A )

20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
