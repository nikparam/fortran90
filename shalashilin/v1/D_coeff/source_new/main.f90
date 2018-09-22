PROGRAM Dcoeff

	IMPLICIT NONE

	INTEGER :: i, j, astatus, &
		   ITOL, ITASK, ISTATE, IOPT, LZW, LRW, IWORK(30), LIW, MF, IPAR, LRP, &
		   switch, key, npts, NumG

	LOGICAL :: exists_ne, exists_D

	CHARACTER :: fname*100, x1*4, file_params*40

	DOUBLE PRECISION :: t_start, t_finish, &
			    T, TOUT, MAXT, &
			    RTOL, ATOL, &
			    junk(6), params(15), m, Step, &
			    N, E, LE, dE

 	DOUBLE COMPLEX :: DUMMY(1,1)

	DOUBLE PRECISION, ALLOCATABLE :: x(:), wts(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:), q(:), p(:), omega(:), phase(:)

	DOUBLE COMPLEX, ALLOCATABLE :: Y(:), ksi(:), eta(:), &
				       S(:,:), H(:,:), &
				       M1(:,:), L(:,:), &
				       ZWORK(:), RPAR(:), &
				       R(:,:) 

	EXTERNAL :: FEXD

	OPEN(UNIT = 10, FILE = 'init_cond.txt', STATUS = 'OLD',&
	     FORM = 'FORMATTED', ACTION = 'READ')
	READ(10,*) switch, key
	READ(10,*) NumG, MAXT, Step
	READ(10,*) dE
	READ(10,*) m

	ALLOCATE(Y(NumG), omega(NumG), phase(NumG), STAT=astatus)
	ALLOCATE(S(NumG,NumG), H(NumG,NumG), M1(NumG,NumG), L(NumG,NumG), R(NumG,NumG), STAT=astatus)

	DO i=1,NumG
		READ(10,*) omega(i), Y(i)
	END DO	
	CLOSE(UNIT = 10)

	file_params = '../potential/params.txt'
	CALL read_params(switch,file_params,params)

	T = 0.0D0
	TOUT = Step
	ITOL = 1
	RTOL = 1.0D-13
	ATOL = 1.0D-13
	ITASK = 1
	ISTATE = 1
	IOPT = 0
	MF = 20
	LZW = 8 * NumG + 2 * NumG**2
	LRW = 20 + NumG
	LIW = 30 + NumG
	LRP = NumG**2

	ALLOCATE(ZWORK(LZW), RWORK(LRW), RPAR(LRP), STAT=astatus)

	ALLOCATE(q(NumG), p(NumG), ksi(NumG), eta(NumG), STAT=astatus)
	OPEN(UNIT=10, FILE = '../mpi/init_cond.txt', &
	     STATUS = 'OLD', FORM = 'FORMATTED', ACTION = 'READ')
	READ(10,*) junk(1)
	READ(10,*) junk(1:3)
	READ(10,*) junk(1)
	DO i=1, NumG
		READ(10,*) q(i), p(i)
	END DO
	CLOSE(UNIT = 10)

	phase(1:NumG) = 0.25D0 * DLOG(omega(1:NumG) / (4.0D0 * DATAN(1.0D0)))
	CALL change_var( NumG, q, p, ksi, eta, omega, phase )
	CALL normalize( NumG, ksi, eta, omega, Y )

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

	CALL CPU_TIME(t_start)

	SELECT CASE ( key )
		CASE(0)
			npts = 1
			ALLOCATE( x(1:npts), wts(1:npts) )
			x(:) = 0.0
			wts(:) = 0.0
		CASE(1)
			npts = 100
			ALLOCATE( x(1:npts), wts(1:npts) )
			CALL p_quadrature_rule ( npts, x, wts )	
!		CASE(2)
!			npts = 100
!			ALLOCATE( x(npts), wts(npts) )
!			CALL h_quadrature_rule ( npts, x, wts )	
	END SELECT

	DO 45 WHILE ( TOUT .LE. MAXT )

		DO j = 1, NumG
			WRITE(x1, '(I2.2)') j
			fname = '../mpi/out' // TRIM(x1) // '.out'
			OPEN(UNIT = 10*j, FILE = fname, STATUS = 'OLD', &
		             FORM = 'FORMATTED', ACTION = 'READ')
			READ(10*j,*) junk(1), q(j), p(j)
		END DO

		CALL change_var( NumG, q, p, ksi, eta, omega, phase )
		CALL overlap(NumG, ksi, eta, omega, S)
		CALL hamiltonian(NumG, key, ksi, eta, omega, params, npts, x, wts, dE, S, H, L, M1)
		CALL coeff_matrix( NumG, m, params, omega, ksi, S, H, M1, R )

		N = DOT_PRODUCT( Y, MATMUL( S, Y ) )
		E = DOT_PRODUCT( Y, MATMUL( H, Y ) ) / N
		LE = DOT_PRODUCT( Y, MATMUL( L, Y ) ) / N
		WRITE(15,35) T, N, E, LE

		RPAR(1:LRP) = RESHAPE( R, (/ NumG**2 /) )
		WRITE(16,25) T, Y(:), (DBLE(Y(:))**2 + DIMAG(Y(:))**2), SUM(DBLE(Y(:))**2 + DIMAG(Y(:))**2)

		ZWORK(:) = (0.0D0, 0.0D0)
		RWORK(:) = 0.0D0
		IWORK(:) = 0
		ISTATE = 1

		CALL ZVODE(FEXD, NumG, Y, T, TOUT, ITOL, RTOL, &
			   ATOL, ITASK, ISTATE, IOPT, &
			   ZWORK, LZW, RWORK, LRW, IWORK, LIW, &
			   DUMMY, MF, RPAR, IPAR)

		S = ( 0.0D0, 0.0D0 )
		H = ( 0.0D0, 0.0D0 )
		M1 = ( 0.0D0, 0.0D0 )
		L = ( 0.0D0, 0.0D0 )

25		FORMAT(F16.6,1X,100F16.6)
35		FORMAT(F16.6,1X,F16.8,'  ',F16.8,' ',F16.8)
		IF (ISTATE .LT. 0) GOTO 85

45	TOUT = TOUT + Step

	CALL CPU_TIME(t_finish)

	CLOSE( UNIT = 15 )
	CLOSE( UNIT = 16 )

	WRITE(6,55) ( t_finish - t_start )
55	FORMAT('CPU time= ', F10.4)


		STOP
85		WRITE(6,95) ISTATE
95		FORMAT(///'  Error halt: ISTATE =',I3)
END PROGRAM Dcoeff
