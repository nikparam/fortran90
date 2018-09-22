SUBROUTINE coeff_matrix( NumG, m, params, omega, ksi, coeff, E, dE, S, H, H1, H2, M1, R )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: m, params(15), omega(NumG), E, dE
	DOUBLE COMPLEX, INTENT(IN) :: ksi(NumG), coeff(NumG), S(NumG,NumG), &
				      H(NumG,NumG), H1(NumG,NumG), H2(NumG,NumG), M1(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: R(NumG,NumG)

	INTEGER :: i, j
	DOUBLE PRECISION :: T, dV(NumG)
	DOUBLE COMPLEX :: dksi(NumG), deta(NumG), Yvec(NumG), &
			  TEMP1(NumG,NumG), TEMP2(NumG,NumG), TEMP3, corr(NumG,NumG), &
			  SI(NumG,NumG), &
			  z_zdot(NumG,NumG)

	DO i=1,NumG
		CALL diff_potential_energy( DBLE( ksi(i) ) / omega(i) , params, dV(i) )
	END DO

	CALL inverse(NumG, S, SI)

	TEMP1 = MATMUL( H, MATMUL( SI, M1 ) )
	TEMP2 = CONJG( TRANSPOSE( TEMP1 ) )
	TEMP3 = DOT_PRODUCT( coeff, MATMUL( H, MATMUL( SI, MATMUL( H, coeff ) ) ) )


	Yvec(:) = ( 0.0D0, 0.0D0 )

	DO i = 1, NumG
		DO j = 1, NumG
			Yvec(i) = Yvec(i) + &
				  CONJG( coeff(j) ) * coeff(i) * ( H1(j,i) - TEMP1(j,i) ) + &
				  coeff(j) * CONJG( coeff(i) ) * ( H2(i,j) - TEMP2(i,j) )
		END DO
	END DO

	dksi(1:NumG) = (/ ( DCMPLX( omega(i) * DIMAG( ksi(i) ) / m, -dV(i) ), i=1,NumG ) /)
	deta(1:NumG) = (/ ( DCMPLX( -DBLE( ksi(i) ) * DIMAG( ksi(i) ), 0.0 ), i=1,NumG ) /)

	z_zdot(:,:) = M1 * SPREAD( dksi, 1, NumG ) + S * SPREAD( deta, 1, NumG )

	corr(:,:) = 2 * DBLE( DOT_PRODUCT( Yvec, dksi ) ) * ( H - E * S ) / ( TEMP3 - E * E )

	R = -( 0.0D0, 1.0D0 ) *  MATMUL( SI, H - (0.0D0, 1.0D0) * z_zdot  - (0.0D0, 1.0D0) * corr )

20 	FORMAT( 100("("F14.6,SP,F14.6,"i) ") )

	RETURN
END SUBROUTINE coeff_matrix
