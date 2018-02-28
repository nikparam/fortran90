PROGRAM main

	IMPLICIT NONE
	INTEGER :: switch, i, NumG, info, astatus
	DOUBLE PRECISION :: junk(6), params(15), Step, m, MAXT
	DOUBLE PRECISION, ALLOCATABLE :: q(:), p(:)
	CHARACTER (LEN = 40) :: file_params

	OPEN(UNIT = 10, FILE = 'init_cond.txt', STATUS = 'OLD',&
	     FORM = 'FORMATTED', ACTION = 'READ')
	READ(10,*) switch
	READ(10,*) NumG, MAXT, Step
	READ(10,*) m
	ALLOCATE( q(NumG), p(NumG), STAT = astatus )
	DO i=1, NumG
		READ(10,*) q(i), p(i)
	END DO
	
 	file_params = '../potential/params.txt'
	CALL read_params(switch, file_params, params)

!	WRITE(*,*) NumG, NSteps, m
!	DO i=1, NumG
!		WRITE(6,20) q(i), p(i)
!	END DO
!20	FORMAT(F14.8,1X,F14.8)

	CALL trajectories(NumG,q(1:NumG),p(1:NumG),MAXT,Step,m,params,info)
	IF (info .EQ. 1) THEN
		WRITE(*,*) '!!SUCCESS!!'
	END IF
END PROGRAM main

