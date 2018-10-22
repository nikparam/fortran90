PROGRAM Dcoeff

	IMPLICIT NONE

	INTEGER :: i, j, astatus, &
		   ITOL, ITASK, ISTATE, IOPT, LZW, LRW, LIW, MF, IPAR(4), LRP, &
		   switch, key, npts, NumG, NEQ, mean

	LOGICAL :: exists_ne, exists_D, exists_qp, exists_coords, exists_w

	CHARACTER :: fname*100, x1*4, file_params*40

	DOUBLE PRECISION :: t_start, t_finish, &
			    T, TOUT, MAXT, &
			    RTOL, ATOL, &
			    junk(6), params(15), m, Step, V, &
			    N, E, LE, KE, VE, &
			    dE, q_0, q_0_square, p_0, width, lambda

 	DOUBLE COMPLEX :: DUMMY(1,1)

	INTEGER, ALLOCATABLE :: IWORK(:)

	DOUBLE PRECISION, ALLOCATABLE :: x(:), wts(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:), q(:), p(:), omega(:), phase(:), RPAR(:)

	DOUBLE COMPLEX, ALLOCATABLE :: Y(:), ksi(:), eta(:), D(:),&
				       S(:,:), H(:,:), &
				       M1(:,:), M2(:,:), L(:,:), &
				       ZWORK(:), &
				       R(:,:), tmp(:) 

	EXTERNAL :: FEXALL

	OPEN(UNIT = 10, FILE = 'init_cond.txt', STATUS = 'OLD',&
	     FORM = 'FORMATTED', ACTION = 'READ')
	READ(10,*) mean
	READ(10,*) switch, key
	READ(10,*) lambda
	READ(10,*) NumG, MAXT, Step
	READ(10,*) dE
	READ(10,*) m

	ALLOCATE(Y(3 * NumG), omega(NumG), phase(NumG), STAT=astatus)
	ALLOCATE(q(NumG), p(NumG), D(NumG), ksi(NumG), eta(NumG), STAT=astatus)
	ALLOCATE(S(NumG,NumG), H(NumG,NumG), M1(NumG,NumG), M2(NumG,NumG), L(NumG,NumG), R(NumG,NumG), STAT=astatus)
	ALLOCATE(tmp(NumG))


	DO i=1,NumG
		READ(10,*) omega(i), Y(2*NumG+i), q(i), p(i)
	END DO	
	CLOSE(UNIT = 10)

	file_params = '../potential/params.txt'
	CALL read_params(switch,file_params,params)

	T = 0.0D0
	TOUT = Step
	ITOL = 1
	RTOL = 1.0D-14
	ATOL = 1.0D-14
	ITASK = 1
	ISTATE = 1
	IOPT = 1
	MF = 20
	NEQ = 3 * NumG
	LZW = 8 * NEQ + 2 * NEQ**2
	LRW = 20 + NEQ
	LIW = 30 + NEQ

	ALLOCATE(ZWORK(LZW), RWORK(LRW), IWORK(LIW), STAT=astatus)

	phase(1:NumG) = 0.25D0 * DLOG(m * omega(1:NumG) / (4.0D0 * DATAN(1.0D0)))

	IF ( mean .GT. 0.0D0 .OR. dE .GT. 0.D0 ) THEN
		p(1:NumG) =  DSQRT( 2 * m * dE )
	END IF

	CALL change_var( NumG, q, p, ksi, eta, m, omega, phase )
	CALL normalize( NumG, ksi, eta, m, omega, Y(2*NumG+1:3*NumG) )

	INQUIRE(FILE = "qp_mean.out", EXIST = exists_qp)
	IF ( exists_qp ) THEN
		OPEN(UNIT = 14, FILE = 'qp_mean.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 14, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 14, FILE = 'qp_mean.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	INQUIRE(FILE = "norm_energy.out", EXIST = exists_ne)
	IF ( exists_ne ) THEN
		OPEN(UNIT = 15, FILE = 'norm_energy.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 15, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 15, FILE = 'norm_energy.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	INQUIRE(FILE = "D_coeff.out", EXIST = exists_D)
	IF ( exists_D ) THEN
		OPEN(UNIT = 16, FILE = 'D_coeff.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 16, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 16, FILE = 'D_coeff.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	INQUIRE(FILE = "coords.out", EXIST = exists_coords)
	IF ( exists_coords ) THEN
		OPEN(UNIT = 17, FILE = 'coords.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 17, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 17, FILE = 'coords.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	INQUIRE(FILE = "width.out", EXIST = exists_w)
	IF ( exists_w ) THEN
		OPEN(UNIT = 18, FILE = 'width.out', FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 18, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 18, FILE = 'width.out', FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	CALL CPU_TIME(t_start)

	npts = 100
	ALLOCATE( x(1:npts), wts(1:npts) )
	CALL p_quadrature_rule ( npts, x, wts )	
	
	Y(1:NumG) = (/ ( DCMPLX( q(i), 0.0D0 ), i = 1, NumG ) /)
	Y(NumG+1:2*NumG) = (/ ( DCMPLX( p(i), 0.0D0 ), i = 1, NumG ) /)
	IPAR(1) = NumG
	IPAR(2) = key
	IPAR(3) = mean
	IPAR(4) = npts

	LRP = 2*(NumG+npts)+16
	ALLOCATE( RPAR(LRP) )
	RPAR(1:NumG) = omega(1:NumG)
	RPAR(NumG+1:2*NumG) = phase(1:NumG)
	RPAR(2*NumG+1) = m
	RPAR(2*NumG+2:2*NumG+16) = params(1:15)
	RPAR(2*NumG+17:2*NumG+16+npts) = x(1:npts)
	RPAR(2*NumG+17+npts:2*NumG+16+2*npts) = wts(1:npts)
	RPAR(2*NumG+17+2*npts) = lambda

	IWORK(6) = 1000
	DO 55 WHILE ( TOUT .LE. MAXT )

		ISTATE = 1

		CALL ZVODE(FEXALL, 3*NumG, Y, T, TOUT, ITOL, RTOL, &
			   ATOL, ITASK, ISTATE, IOPT, &
			   ZWORK, LZW, RWORK, LRW, IWORK, LIW, &
			   DUMMY, MF, RPAR, IPAR)

		q(1:NumG) = DBLE( Y(1:NumG) )
		p(1:NumG) = DBLE( Y(NumG+1:2*NumG) )
		D(1:NumG) = Y(2*NumG+1:3*NumG)

		CALL change_var( NumG, q, p, ksi, eta, m, omega, phase )
		CALL overlap( NumG, ksi, eta, m, omega, S )
		CALL hamiltonian( NumG, key, ksi, eta, m, omega, params, npts, x, wts, S, H, L, M1, M2 )

		N = DOT_PRODUCT( D, MATMUL( S, D ) )
		E = DOT_PRODUCT( D, MATMUL( H, D ) ) / N
		LE = DOT_PRODUCT( D, MATMUL( L, D ) ) / N
		KE = 0.5 * DOT_PRODUCT( D, MATMUL( H + L, D ) ) / N
		VE = 0.5 * DOT_PRODUCT( D, MATMUL( H - L, D ) ) / N
		q_0 = DOT_PRODUCT( D, MATMUL( M1, D ) )
		CALL potential_energy(q_0,params,V)
		p_0 = DOT_PRODUCT( D, MATMUL( -DCMPLX(0.0D0,1.0D0) * ( - m * SPREAD( omega, 1, NumG ) * M1 + &
									     SPREAD( ksi, 1, NumG ) * S ), D ) )
		q_0_square = DOT_PRODUCT( D, MATMUL( M2, D ) )
		width = q_0_square - q_0**2
		WRITE(14,'(3F16.6)') T, q_0, p_0

		tmp = MATMUL( S, D )

		WRITE(15,25) T, &
			     N, &
			     E, &
			     KE, &
			     VE, &
			     LE, &
			     SUM( (/ ( DBLE( H(i,i) ), i=1,NumG ) /) )
		WRITE(17,35) T, q, p
		WRITE(16,45) T, D, ABS(D)
		WRITE(18,'(2F16.6)') T, width

25		FORMAT(F16.6,1X,7F16.8)
35 		FORMAT(F16.6,10F16.6,'      ',10F16.6)
45		FORMAT(F16.6,1X,100F16.6)
		IF (ISTATE .EQ. -1) THEN
			ISTATE = 2
			IWORK(6) = 2 * IWORK(6)
		END IF
		IF (ISTATE .LT. -1) GOTO 85

55	TOUT = TOUT + Step

	CALL CPU_TIME(t_finish)

	CLOSE( UNIT = 15 )
	CLOSE( UNIT = 16 )
	CLOSE( UNIT = 17 )
	CLOSE( UNIT = 18 )

	WRITE(6,65) ( t_finish - t_start )
65	FORMAT('CPU time= ', F10.4)


		STOP
85		WRITE(6,95) ISTATE
95		FORMAT(///'  Error halt: ISTATE =',I3)
END PROGRAM Dcoeff
