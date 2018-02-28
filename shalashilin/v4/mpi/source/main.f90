PROGRAM main

	INCLUDE 'mpif.h'

	INTEGER :: switch, i, NumG, info, astatus
	DOUBLE PRECISION :: junk(6), params(15), Step, m, MAXT, ref_num
	DOUBLE PRECISION, ALLOCATABLE :: q(:), p(:)
	CHARACTER (LEN = 40) :: file_params

	INTEGER, PARAMETER :: send_data_tag = 2001, return_data_tag = 2002

	INTEGER :: my_id, root_process, ierr, STATUS(MPI_STATUS_SIZE)
	INTEGER :: num_procs, an_id, num_rows_to_receive
	INTEGER :: avg_rows_per_process, num_rows_to_send
	INTEGER :: start_row, end_row

	DOUBLE PRECISION, ALLOCATABLE :: q2(:), p2(:)

	root_process = 0

	CALL MPI_INIT(ierr)

	CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr)
	CALL MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr)

	if ( my_id .EQ. root_process ) THEN

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

		avg_rows_per_process = NumG / num_procs

		DO an_id = 1, num_procs - 1

			start_row = an_id * avg_rows_per_process + 1
			end_row = start_row + avg_rows_per_process - 1

			IF ( an_id .EQ. ( num_procs - 1 ) ) THEN
				end_row = NumG
			END IF

			num_rows_to_send = end_row - start_row + 1

			CALL MPI_SEND( num_rows_to_send, 1, MPI_INT, &
				       an_id, send_data_tag, MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( avg_rows_per_process, 1, MPI_INT, &
				       an_id, send_data_tag, MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( q(start_row), num_rows_to_send, &
				       MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( p(start_row), num_rows_to_send, &
				       MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( MAXT, 1, MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( Step, 1, MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( m, 1, MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( params, 15, MPI_DOUBLE_PRECISION, an_id, &
				       send_data_tag, MPI_COMM_WORLD, ierr)

		END DO

		CALL trajectories(avg_rows_per_process, avg_rows_per_process, &
				  my_id, &
				  q(1:avg_rows_per_process), &
				  p(1:avg_rows_per_process), &
				  MAXT, Step, m, params)
	ELSE

		CALL MPI_RECV( num_rows_to_receive, 1, MPI_INT, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( avg_rows_per_process, 1, MPI_INT, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		ALLOCATE( q2(num_rows_to_receive), p2(num_rows_to_receive) )

		CALL MPI_RECV( q2, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			       root_process, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( p2, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			       root_process, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( MAXT, 1, MPI_DOUBLE_PRECISION, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( Step, 1, MPI_DOUBLE_PRECISION, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( m, 1, MPI_DOUBLE_PRECISION, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL MPI_RECV( params, 15, MPI_DOUBLE_PRECISION, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		CALL trajectories(num_rows_to_receive, avg_rows_per_process, &
				  my_id, &
				  q2, p2, &
				  MAXT, Step, m, params)

	END IF

	CALL MPI_FINALIZE(ierr)

END PROGRAM main

