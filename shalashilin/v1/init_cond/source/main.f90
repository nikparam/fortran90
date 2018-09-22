PROGRAM main

	USE GRID

	INTEGER :: switch, NumV, NumG, astatus
	DOUBLE PRECISION :: params(15), omega, m, factor, q_lim, TMAX, TSTEP, V_junk, N, E
	DOUBLE PRECISION, ALLOCATABLE :: V_lim(:), q_array(:), p_array(:), omega_array(:), phase_array(:)
	DOUBLE COMPLEX, ALLOCATABLE :: eigen_states(:), eigen_vectors(:,:)
	CHARACTER :: file_params*40
	CHARACTER :: qpname*40, ename*40, sname*40
	LOGICAL :: exists_qp, exists_e, exists_s

	OPEN(UNIT = 10, FILE = 'init_cond.txt', STATUS = 'OLD', &
	     FORM = 'FORMATTED', ACTION = 'READ')

	READ(10,*) switch
	READ(10,*) TMAX, TSTEP
	READ(10,*) omega
	READ(10,*) m
	READ(10,*) factor
	READ(10,*) q_lim
	READ(10,*) NumV

	ALLOCATE(V_lim(NumV), STAT = astatus)

	DO i = 1, NumV
		READ(10,*) V_lim(i)
	END DO

	CLOSE(UNIT = 10)

	file_params = '../potential/params.txt'
	CALL read_params(switch, file_params, params)

	CALL get_grid( omega, m, factor, params, q_lim, V_lim, NumV, q_array, p_array, omega_array, NumG)

	qpname = 'qp_init.out'
	INQUIRE(FILE = qpname, EXIST = exists_qp)
	IF ( exists_qp ) THEN
		OPEN(UNIT = 20, FILE = qpname, FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 20, STATUS = 'DELETE')
	END IF
	OPEN(UNIT = 20, FILE = qpname, FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	WRITE(20,'(I2.2)') switch
	WRITE(20,45) NumG, TMAX, TSTEP
	WRITE(20,'(F12.6)') m

	DO i = 1, NumG
		CALL potential_energy( q_array(i), params, V_junk )
		WRITE(20,35) q_array(i), p_array(i), 0.5D0 * p_array(i)**2 / m + V_junk
	END DO

	CLOSE(UNIT = 20)

35	FORMAT(F14.8,' ', F14.8,' ', F14.8,' ',F14.8)
45	FORMAT(I4.2,' ', F14.8,' ', F14.8)

	ALLOCATE(eigen_states(NumG), eigen_vectors(NumG, NumG), phase_array(NumG), STAT = astatus)
	phase_array(1:NumG) = 0.25D0 * DLOG(omega_array(1:NumG) / (4.0D0 * DATAN(1.0D0)))
	CALL GET_STATES(NumG, omega_array, phase_array, params, q_array, p_array, eigen_states, eigen_vectors)
	CALL sort( eigen_states, NumG, eigen_vectors )

	sname = 'eigen_vectors.out'
	INQUIRE(FILE = sname, EXIST = exists_s)
	IF ( exists_s ) THEN
		OPEN(UNIT = 20, FILE = sname, FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 20, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 20, FILE = sname, FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	WRITE(20,'(I2.2)') NumG

	WRITE(20,'(100F12.8)') DBLE(eigen_states)

	DO i = 1, NumG
		WRITE(20,55) eigen_vectors(i,:)
	END DO

55	FORMAT( 100( F19.8,' ', F19.8 ) )

	CLOSE(UNIT = 20)

	sname = 'format_output.out'
	INQUIRE(FILE = sname, EXIST = exists_s)
	IF ( exists_s ) THEN
		OPEN(UNIT = 20, FILE = sname, FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 20, STATUS = 'DELETE')
	END IF

	OPEN(UNIT = 20, FILE = sname, FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	WRITE(20,*) "	q		p		E"
	DO i = 1, NumG
		CALL potential_energy( q_array(i), params, V_junk )
		WRITE(20,'(F14.6," ",F14.6," ",F14.6)') q_array(i), p_array(i), 0.5D0 * p_array(i)**2 / m + V_junk
	END DO

	WRITE(20,*)

	DO i = 1, NumG
		CALL norm_energy( NumG, q_array, p_array, eigen_vectors(:,i), N, E, omega_array, params, phase_array )
		WRITE(20,*) "------------------------------------------------------------------------"
		WRITE(20,'("Energy[ ",I2.2," ]=",F14.6," and ",F14.6", Norm = ",F14.6)') i, DBLE( eigen_states(i) ), E, N
		WRITE(20,*) "-------------------------------------------------------------------------"
		DO j = 1, NumG
			WRITE(20,'("D[",I2.2,",",I2.2,"]=   (",F14.6,",",SP,F14.6,")")') j, i, eigen_vectors(j,i)
		END DO
		WRITE(20,*)
	END DO

	ename = 'energy_states.out'
	INQUIRE(FILE = ename, EXIST = exists_e)
	IF ( exists_e ) THEN
		OPEN(UNIT = 20, FILE = ename, FORM = 'FORMATTED', &
		     STATUS = 'OLD', ACTION = 'WRITE')
		CLOSE(UNIT = 20, STATUS = 'DELETE')
	END IF
	OPEN(UNIT = 20, FILE = ename, FORM = 'FORMATTED', &
	     STATUS = 'NEW', POSITION = 'APPEND', ACTION = 'WRITE')

	DO i = 1, NumG
		WRITE(20,'(4F16.8)') omega * ( i - 0.5D0 ), &
				     DBLE( eigen_states(i) ), &
				     ABS( omega * ( i - 0.5D0 ) - DBLE( eigen_states(i) ) ), &
				     1.0D2 * (ABS( omega * ( i - 0.5D0 ) - DBLE( eigen_states(i) ) ) ) / omega / (i - 0.5D0)
	END DO

	CLOSE(UNIT = 20)

	CONTAINS
		SUBROUTINE sort(array1, size_array, array2 )

			INTEGER, INTENT(IN) :: size_array
			DOUBLE COMPLEX, INTENT(INOUT) :: array1( size_array ), array2( size_array, size_array )
			DOUBLE COMPLEX :: vector( size_array )
			INTEGER :: i, j, k
			DOUBLE COMPLEX :: element

			DO i = 1, size_array
				DO j = i, size_array
					IF ( DBLE( array1(j) ) .LT. DBLE( array1(i) ) ) THEN
						element = array1(j)
						array1(j) = array1(i)
						array1(i) = element
						DO k = 1, size_array
							vector(k) = array2(k,j)
							array2(k,j) = array2(k,i)
							array2(k,i) = vector(k)
						END DO
					END IF
				END DO
			END DO

			RETURN

		END SUBROUTINE sort

END PROGRAM main
