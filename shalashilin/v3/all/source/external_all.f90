SUBROUTINE FEXALL( N, T, Y, YDOT, RPAR, IPAR )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: N, IPAR(4)
	DOUBLE PRECISION, INTENT(IN) :: T, RPAR(2*IPAR(1)+17+2*IPAR(4))
	DOUBLE COMPLEX, INTENT(IN) :: Y(N)
	DOUBLE COMPLEX, INTENT(OUT) :: YDOT(N)

	INTEGER :: i, j, mean, key, NumG, npts, IPIV( ipar(1) ), INFO
	DOUBLE PRECISION :: m, params(15), lambda, &
			    Norm_D

	DOUBLE PRECISION :: omega( IPAR(1) ), phase( IPAR(1) ), &
			    x( IPAR(4) ), wts( IPAR(4) ), &
			    q( IPAR(1) ), p( IPAR(1) ), &
			    dV( IPAR(1) ), dV_mean

	DOUBLE COMPLEX :: ksi( IPAR(1) ), eta( IPAR(1) ), &
			  S( IPAR(1), IPAR(1) ), H( IPAR(1), IPAR(1) ), &
			  L( IPAR(1), IPAR(1) ), tmp, &
			  M1( IPAR(1), IPAR(1) ), M2( IPAR(1), IPAR(1) ), &
			  R( IPAR(1), IPAR(1) ), D( IPAR(1) ), dV_matrix(ipar(1),ipar(1))

	INTERFACE
		SUBROUTINE potential_energy( x, params, V )
			DOUBLE PRECISION, INTENT(IN) :: x, params(15)
			DOUBLE PRECISION, INTENT(OUT) :: V
		END SUBROUTINE potential_energy
		SUBROUTINE diff_potential_energy( x, params, V )
			DOUBLE PRECISION, INTENT(IN) :: x, params(15)
			DOUBLE PRECISION, INTENT(OUT) :: V
		END SUBROUTINE diff_potential_energy

	END INTERFACE

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

	IF ( mean .EQ. 0 ) THEN
		CALL quadrature( NumG, diff_potential_energy, ksi, eta, x, wts, npts, m, omega, params, dV_matrix )
		DO i = 1, NumG
			dV(i) = DBLE( dV_matrix(i,i) )
		END DO
	ELSE IF( mean .EQ. 1 ) THEN
		CALL quadrature( NumG, diff_potential_energy, ksi, eta, x, wts, npts, m, omega, params, dV_matrix )
		dV_mean = DOT_PRODUCT( D, MATMUL( dV_matrix, D ) )
		dV = (/ ( dV_mean, i=1,NumG ) /)
	ELSE IF ( mean .EQ. 2 ) THEN
		CALL quadrature( NumG, diff_potential_energy, ksi, eta, x, wts, npts, m, omega, params, dV_matrix )
		dV_mean = SUM( (/ ( dV_matrix(i,i), i = 1, NumG ) /) )
		dV = (/ ( dV_mean / NumG, i = 1, NumG ) /)
	ELSE
		CALL quadrature( NumG, diff_potential_energy, ksi, eta, x, wts, npts, m, omega, params, dV_matrix )
		Norm_D = SUM( (/ ( CONJG(D(i))*D(i), i=1,NumG ) /) )
		dV_mean = SUM( (/ ( dV_matrix(i,i) * CONJG(D(i)) * D(i), i=1,NumG ) /) )
		dV = (/ ( dV_mean / Norm_D, i=1,NumG ) /)
	
	END IF

	CALL coeff_matrix( NumG, m, params, omega, lambda, ksi, dV, S, H, M1, R )

	YDOT(1:NumG) = p(1:NumG) / m
	YDOT(NumG+1:2*NumG) = -dV(1:NumG)
	YDOT(2*NumG+1:3*NumG) = MATMUL( R, D )

END SUBROUTINE FEXALL
