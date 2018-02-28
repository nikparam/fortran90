PROGRAM example

!	Example for DLSODE integrator from ODEPACK
!
!	Compile with command:
!	
!	gfortran -o example.o example.f90 external.f90 opkdmain.f opkda1.f opkda2.f
!
!	where example.f90 -- current script
!	      external.f90 -- script with subroutines for derivative vector and Jacobian matrix
!	      opkd -- ODEPACK library files
!
!	while compiling, there will be certain warnings. Discard them without second thought

	
	IMPLICIT NONE
	EXTERNAL ::  FEX, JEX ! names of subroutines in external file

!	FEX = vector of derivatives
!
!	JEX = Jacobian matrix

	INTEGER :: NEQ, IOUT, ITOL, ITASK, ISTATE, IOPT, LRW, MF, IWORK(23), LIW ! INTEGER parameters

!	NEQ = number of equation
!
!	IOUT = number of steps
!
!	ITOL = 1 if ATOL is scalar, 2 if ATOL is array
!
!	ITASK = flag, indicating task of DLSODE, 1 for normal calculation of Y ar point TOUT
!
!	ISTATE = index to specify a state of calculation:
! 		1 if its a first call
!		2 if subsequent call
!		If ITASK < 0, there is an error --> look in opkdmain.f 
!
!	IOPT = flag for optional input, 0 if none
!
!	LRW = length of RWORK
!		LRW = 20 + 16 * NEQ				for MF = 10;
!		LRW = 22 + 9 * NEQ + NEQ**2			for MF = 21 or 22;
!		LRW = 22 + 10 * NEQ + ( 2*ML + MU ) * NEQ	for MF = 24 or 25.
!
!	MF = method flag;
!		(10)	Nonstiff (Adams) method, no Jacobian used
!		(21)	Stiff (BDF) method, user-supplied full Jacobian
!		(22)	Stiff method, internally generated full Jacobian
!		(24)	Stiff method, user-supplied banded Jacobian
!		(25) 	Stiff method, internally generated banded Jacobian
!
!	IWORK(LIW) = integer work array;
!
!	LIW = length of IWORK
!		LIW = 20 for MF = 20
!		LIW = 20 + NEQ for MF = 21, 22, 24 or 25
 
	DOUBLE PRECISION :: ATOL(3), RWORK(58), RTOL, T, TOUT, Y(3) ! REAL*8 parameters
!	ATOL = absolute tolerance, scalar or array
!
!	RWORK = real work array
!
!	RTOL = relative tolerance
!
!	T = value of variable; on return it will be current value of T (normally, TOUT)
!
!	TOUT = next point
!
!	Y = array of values of function at time step T


! Number of equations to solve:
	NEQ = 3

! Initial conditions:
	Y(1) = 1.D0
	Y(2) = 0.D0
	Y(3) = 0.D0
	T = 0.D0

! Next point at time:
	TOUT = 0.4D0

! ATOL -- array, RTOL -- scalar, tolerance will depend only upon absolute tolerance:
	ITOL = 2
	RTOL = 1.D-4
	ATOL(1) = 1.D-6
	ATOL(2) = 1.D-10
	ATOL(3) = 1.D-6

! Setting parameters of calculation:
	ITASK  = 1
	ISTATE = 1
	IOPT = 0
	LRW = 58
	LIW = 23
	MF = 21

! Solve system of ODE from T for 12 values of TOUT, starting from TOUT = 0.4, each time multiplying it by 10:

	DO 40 IOUT = 1,12 ! Nice F77 way to iterate: do condition ( label 40 ) while IOUT in [1,12, step = 1]

!		Call DLSODE subroutine from opkdmain.f

		CALL DLSODE(FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, ISTATE, &
			   IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)

!		Write current time and values of function at current time

		WRITE(6,20) T, Y(1), Y(2), Y(3) ! 6 -- output with WRITE to the screen, 20 -- label, point on format
		20 FORMAT (' AT T =',D12.4,'  Y =',3D14.6) ! This will write one R*8 number with 4 digits after point and 12 digits total
							   !	      and three R*8 numbers with 6 digits after point and 14 digits total

!		Check, whether there are any errors;
!			If yes ( ISTATE < 0 ), write them (See 80)

		IF ( ISTATE .LT. 0 ) GO TO 80

	40 TOUT = TOUT * 10.D0 ! condition, that will be executed, while IOUT iterates through [1,12,step = 1]

! write number of steps ( IWORK(11) ),
!	number of derivative evaluations ( IWORK(12) ),
!	number of Jacobian evaluations ( IWORK(13) ).

	WRITE(6,60) IWORK(11), IWORK(12), IWORK(13)
	60 FORMAT(/' NO. STEPS =',I4',  No. F-S =',I4',  NO. J-S =',I4) ! This will write three integers with total 4 digits
	STOP

	80 WRITE(6,90) ISTATE
	90 FORMAT(///' ERROR HALT.. ISTATE=',I3) ! This will write one integer with 3 digits total
	STOP

END PROGRAM EXAMPLE
