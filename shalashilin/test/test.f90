PROGRAM sum_vec

	INCLUDE 'mpif.h'

	INTEGER, PARAMETER :: send_data_tag = 2001, return_data_tag = 2002

	INTEGER :: my_id, root_process, ierr, STATUS(MPI_STATUS_SIZE)
	INTEGER :: num_procs, an_id, num_rows_to_receive
	INTEGER :: avg_rows_per_process, num_rows, num_rows_to_send
	INTEGER :: start_row, end_row

	DOUBLE PRECISION, ALLOCATABLE :: vector(:), vector2(:)
	DOUBLE PRECISION :: partial_sum, total_sum

	root_process = 0

	CALL MPI_INIT(ierr)

	CALL MPI_COMM_RANK( MPI_COMM_WORLD, my_id, ierr)
	CALL MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr)

	IF ( my_id .EQ. root_process ) THEN

		OPEN(UNIT = 10, FILE = 'data.txt', STATUS = 'OLD', &
		     FORM = 'FORMATTED', ACTION = 'READ')

		READ(10,*) num_rows
		ALLOCATE( vector(num_rows) )
		avg_rows_per_process = num_rows / num_procs

		DO i = 1, num_rows
			READ(10,*) vector(i)
		END DO

		DO an_id = 1, num_procs - 1

			start_row = an_id * avg_rows_per_process + 1
			end_row = start_row + avg_rows_per_process  - 1

			IF ( an_id .EQ. ( num_procs - 1 ) ) THEN
				end_row = num_rows
			END IF

			num_rows_to_send = end_row - start_row + 1

			CALL MPI_SEND( num_rows_to_send, 1, MPI_INT, &
				       an_id, send_data_tag, MPI_COMM_WORLD, ierr)

			CALL MPI_SEND( vector(start_row), num_rows_to_send, &
				       MPI_DOUBLE_PRECISION, an_id, send_data_tag, &
				       MPI_COMM_WORLD, ierr)

		END DO

		total_sum = 0.0D0
		DO i = 1, avg_rows_per_process
			total_sum = total_sum + vector(i)
		END DO

		DO an_id = 1, num_procs - 1

			CALL MPI_RECV( partial_sum, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,&
				       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

			total_sum = total_sum + partial_sum
		END DO

		WRITE(*,*) total_sum

	ELSE

		CALL MPI_RECV( num_rows_to_receive, 1, MPI_INT, root_process, &
			       MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		ALLOCATE(vector2(num_rows_to_receive))

		CALL MPI_RECV( vector2, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			       root_process, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, ierr)

		partial_sum = 0.0D0
		DO i = 1, num_rows_to_receive
			partial_sum = partial_sum + vector2(i)
		END DO

		CALL MPI_SEND( partial_sum, 1, MPI_DOUBLE_PRECISION, root_process, &
			       return_data_tag, MPI_COMM_WORLD, ierr)

	END IF

	CALL MPI_FINALIZE(ierr)


END PROGRAM sum_vec
