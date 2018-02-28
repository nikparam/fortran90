PROGRAM test1

	INCLUDE 'mpif.h'
	INTEGER :: i, ir, is, comm, rank, numproc, ierror, root, request
	INTEGER :: status(MPI_STATUS_SIZE)
	CHARACTER(80) :: message_sent, message_received
	CHARACTER :: x1*3

	message_sent = 'No message sent'
	message_received = 'No message received'

	root = 0

	CALL MPI_INIT(ierror)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierror)

	IF ( numproc .GT. 1 ) THEN
	DO i = 1,5
		IF ( rank .EQ. root ) THEN
				WRITE(x1,'(I2.2)') i
				message_sent =  'Hello from processor 0: ' // x1
				is = i
				CALL MPI_IRECV( ir, 1, MPI_INT, 1, 1, &
					       MPI_COMM_WORLD, request, ierror )

				CALL MPI_SEND( is, 1, MPI_INT, 1, 1, &
					       MPI_COMM_WORLD, ierror )

				CALL MPI_WAIT(request, status, ierror)
		ELSE iF ( rank .EQ. 1 ) THEN
			CALL MPI_RECV( ir, 1, MPI_INT, 0, 1, &
				       MPI_COMM_WORLD, status, ierror )

			message_sent = 'Proc 1 got this message: '//message_received
			is = i * 10
			CALL MPI_SEND( is, 1, MPI_INT, 0, 1, &
				       MPI_COMM_WORLD, ierror)
		END IF
		PRINT *, "Processor ", rank, " sent '",is,"'"
		PRINT *, "Processor ", rank, " received '",ir,"'"
	END DO
	ELSE
		PRINT *, "Not enough processors to demo message passing"
	END IF

	CALL MPI_FINALIZE(ierror)

END PROGRAM test1
