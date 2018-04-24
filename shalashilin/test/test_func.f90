PROGRAM test_func

	INTEGER, ALLOCATABLE :: array(:)
	INTEGER :: N

	N = read_N()

	ALLOCATE( array(N) )

	array(:) = my_func(N)

	WRITE(*,*) array

	CONTAINS

		FUNCTION read_N()

			INTEGER :: read_N

			WRITE(*,'(a)') 'INPUT N: '
			READ(*,*) read_N

		END FUNCTION read_N 

		FUNCTION my_func(N)

				INTEGER :: N
				INTEGER :: my_func(N)

				DO i = 1, N
					my_func(i) = i
				END DO

		END FUNCTION my_func

END PROGRAM test_func
