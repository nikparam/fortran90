SUBROUTINE trajectories(NumG, avg, proc_id, q, p, MAXT, Step, m, params)

	EXTERNAL :: FEXqp
	INTEGER :: avg, i, proc_id, fnum, itraj, NEQ, ITASK, ISTATE, IOPT, LRW, IWORK(50), &
		   LIW, MF, IPAR, ITOL, switch, name1, name2
	DOUBLE PRECISION :: RWORK(100), &
			    DUMMY(1,1), RPAR(16), &
			    ATOL, RTOL, V, dV1, dV2, params(15), m, MAXT
	DOUBLE PRECISION :: q(NumG), p(NumG), Y(2), T, TOUT, Step, X, dV, test
	DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-5, EPS_step = 1.0D-4
	CHARACTER :: fname_traj * 10, fname_en * 20, x1 * 5
	LOGICAL :: exists_traj, exists_en
	REAL :: t_start, t_finish

	NEQ = 2
	ITOL = 1
	RTOL = 1.0D-15
	ATOL = 1.0D-15
	IOPT = 0
	LRW = 100
	LIW = 50
	MF = 20


	RPAR(1) = m
	RPAR(2 : 16) = params(1 : 15)
	IPAR = NumG

	CALL CPU_TIME(t_start)
	DO itraj = 1, NumG
		ISTATE = 1
		ITASK = 1
		T = 0.0D0
		TOUT = Step
		Y(1) = q(itraj)
		Y(2) = p(itraj)
		fnum = itraj + proc_id * avg
		WRITE(x1, '(I2.2)') fnum
		fname_traj = 'out' // TRIM(x1) // '.out'
		fname_en = 'energy' // TRIM(x1) // '.out'
		INQUIRE(FILE = fname_traj, EXIST = exists_traj)
		INQUIRE(FILE = fname_en, EXIST = exists_en)

		name1 = 20 + fnum
		name2 = 30 + fnum

		IF (exists_traj) THEN
			OPEN(UNIT = name1, FILE = fname_traj, FORM = 'FORMATTED', &
			     STATUS = 'OLD', ACTION = 'WRITE')
			CLOSE(UNIT = name1, STATUS = 'DELETE')
		END IF

		IF (exists_en) THEN
			OPEN(UNIT = name2, FILE = fname_en, FORM = 'FORMATTED', &
			     STATUS = 'OLD', ACTION = 'WRITE')
			CLOSE(UNIT = name2, STATUS = 'DELETE')
		END IF

		OPEN(UNIT = name1, FILE = fname_traj, FORM = 'FORMATTED', &
		     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

		OPEN(UNIT = name2, FILE = fname_en, FORM = 'FORMATTED', &
		     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

		DO 13 WHILE ( TOUT .LE. MAXT )

			ISTATE = 1
			RWORK = 0.0D0
			IWORK = 0

			CALL DVODE(FEXqp, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, &
				   IOPT, RWORK, LRW, IWORK, LIW, DUMMY, MF, RPAR, IPAR)

			CALL diff_potential_energy( Y(1), params(1:15), dV )
			IF ( ABS( dV ) .GT. EPS ) THEN
				test = Step * EPS_step / ABS( dV )
			END IF
			WRITE(name1,11) T, Y(:)
			CALL potential_energy( Y(1), params(1:15), V )
			WRITE(name2,11) T, Y(2)**2/(2 * m),  V, Y(2)**2/(2 * m) +  V

11			FORMAT(F12.4,' ',E20.14,' ', E20.14,' ', E20.14)
12			FORMAT(F12.4,:,3F14.6)
			IF ( ISTATE .LT. 0 ) GO TO 15
13		TOUT = TOUT + Step
		CLOSE(UNIT = name1)
		CLOSE(UNIT = name2)
!		WRITE(6,60) IWORK(11), IWORK(12), IWORK(13), IWORK(19), &
!			    IWORK(20), IWORK(21), IWORK(22)
!60		FORMAT(/'  No. steps =',I4,'  No. f-s =',I4, &
!			'  No. J-s =',I4,'  No. LU-s =',I4/ & 
!			'  No. nonlinear iterations =',I4/ &
!			'  No. nonlinear convergence failures =',I4/ &
!			'  No. error test failures =',I4/)
		RWORK = 0.0D0
		IWORK = 0
	END DO		
	CALL CPU_TIME(t_finish)
	WRITE(6,14) t_finish - t_start
14	FORMAT('CPU time= ', F10.4)
!	info = 1
	RETURN

	STOP
15	WRITE(6,16) ISTATE
16	FORMAT(///'  Error halt: ISTATE =',I3)
!	info = -1
	STOP


END SUBROUTINE trajectories
