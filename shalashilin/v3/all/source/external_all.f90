SUBROUTINE FEXALL( N, T, Y, YDOT, RPAR, IPAR )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: N, IPAR(4)
	DOUBLE PRECISION, INTENT(IN) :: T, RPAR(2*IPAR(1)+17+2*IPAR(4))
	DOUBLE COMPLEX, INTENT(IN) :: Y(N)
	DOUBLE COMPLEX, INTENT(OUT) :: YDOT(N)

	INTEGER :: i, mean, key, NumG, npts
	DOUBLE PRECISION :: m, params(15), q_0, lambda, &
			    curvature, xi, xi_prime, rho, z_lambda(ipar(1))

	DOUBLE PRECISION :: omega( IPAR(1) ), phase( IPAR(1) ), &
			    x( IPAR(4) ), wts( IPAR(4) ), &
			    q( IPAR(1) ), p( IPAR(1) ), dV( IPAR(1) )

	DOUBLE COMPLEX :: ksi( IPAR(1) ), eta( IPAR(1) ), &
			  S( IPAR(1), IPAR(1) ), H( IPAR(1), IPAR(1) ), &
			  L( IPAR(1), IPAR(1) ), &
			  M1( IPAR(1), IPAR(1) ), M2( IPAR(1), IPAR(1) ), &
			  R( IPAR(1), IPAR(1) ), D( IPAR(1) ), A( IPAR(1), IPAR(1) ), &
			  tmp( IPAR(1) ), SI( ipar(1), ipar(1) )

	NumG = IPAR(1)
	key = IPAR(2)
	mean = IPAR(3)
	npts = IPAR(4)

	omega(1:NumG) = RPAR(1:NumG)
	phase(1:NumG) = RPAR(NumG+1:2*NumG)
	m = RPAR(2*NumG+1)
	params(1:15) = RPAR(2*NumG+2:2*NumG+16)
	x(1:npts) = RPAR(2*NumG+17:2*NumG+16+npts)
	wts(1:npts) = RPAR(2*NumG+17+npts:2*NumG+16+2*npts)
	lambda = RPAR(2*NumG+17+2*npts)

	q(1:NumG) = (/ ( DBLE( Y(i) ), i = 1, NumG ) /)
	p(1:NumG) = (/ ( DBLE( Y(NumG+i) ), i = 1, NumG ) /)
	D(1:NumG) = (/ ( Y(2*NumG+i), i = 1, NumG ) /)
	
	CALL change_var( NumG, q, p, ksi, eta, m, omega, phase )
	CALL overlap(NumG, ksi, eta, m, omega, S)
	CALL hamiltonian(NumG, key, ksi, eta, m, omega, params, npts, x, wts, S, H, L, M1, M2)

	IF ( mean .EQ. 1 ) THEN
		q_0 = DOT_PRODUCT( D, MATMUL( M1, D ) )
	ELSE IF ( mean .EQ. 2 ) THEN
		q_0 = DOT_PRODUCT( D, MATMUL( S * SPREAD( q, 1, NumG ) , D ) )
	END IF

	IF ( mean .EQ. 0 ) THEN
		DO i = 1, NumG
			CALL diff_potential_energy(q(i), params, dV(i))
		END DO
	ELSE
		DO i = 1, NumG
			CALL diff_potential_energy(q_0, params, dV(i))
		END DO
	END IF	
	CALL coeff_matrix( NumG, mean, m, params, omega, lambda, ksi, q_0, S, H, M1, R )

	YDOT(1:NumG) = p(1:NumG) / m
	YDOT(NumG+1:2*NumG) = -dV(1:NumG)
	YDOT(2*NumG+1:3*NumG) = MATMUL( R, D(1:NumG) )

!	tmp = MATMUL( S, YDOT(2*NumG+1:3*NumG) ) - MATMUL( A, D )
!	WRITE(6,'(3F16.6)') T, DBLE( DOT_PRODUCT( YDOT(2*NumG+1:3*NumG), YDOT(2*NumG+1:3*NumG) ) ), &
!			       DBLE( DOT_PRODUCT( tmp, tmp ) )

!	xi = DBLE( DOT_PRODUCT( YDOT(2*NumG+1:3*NumG), YDOT(2*NumG+1:3*NumG) ) )
!	rho = DBLE( DOT_PRODUCT( tmp, tmp ) )
!	z_lambda = MATMUL( SI, tmp )
!	xi_prime = 4 * DOT_PRODUCT( YDOT(2*NumG+1:3*NumG), z_lambda ) / lambda
!
!	curvature = 2 * xi * rho / xi_prime * ( lambda**2 * xi_prime * rho + &
!						2 * lambda * xi * rho + &
!						lambda**4 * xi * xi_prime ) / ( lambda**2 * xi**2 + rho**2 )**1.5

!	WRITE(6,'(2F20.8)') T, curvature

END SUBROUTINE FEXALL
