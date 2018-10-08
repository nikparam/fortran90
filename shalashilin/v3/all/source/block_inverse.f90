SUBROUTINE block_inverse( N, lambda, M, Mg )

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: N
	DOUBLE PRECISION, INTENT(IN) :: lambda
	DOUBLE COMPLEX, INTENT(IN) :: M(N,N)
	DOUBLE COMPLEX, INTENT(OUT) :: Mg(N,N)

	INTEGER :: LW, RW, IW, INFO, rank, i, j
	DOUBLE PRECISION :: SV(N)
	DOUBLE COMPLEX :: TEMP(N,N), U(N,N), VT(N,N)
	INTEGER, ALLOCATABLE :: IWORK(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
	DOUBLE COMPLEX, ALLOCATABLE :: WORK(:), &
				       A(:,:), B(:,:), C(:,:), &
				       Q(:,:), Ag(:,:), Qg(:,:)


	LW = 3 * N
	RW = 5 * N
	IW = 8 * N

	ALLOCATE( WORK( LW ), RWORK( RW ), IWORK( IW ) )
	TEMP(:,:) = M(:,:)

	CALL ZGESDD( 'N', N, N, TEMP, N, SV, U, N, VT, N, WORK, LW, RWORK, IWORK, INFO )

	rank = 0
	DO i = 1, N
		IF ( SV(i) .GE. lambda ) THEN
			rank = rank + 1
		ELSE
			EXIT
		END IF
	END DO

	IF ( rank .LT. N .AND. rank .GT. 0 ) THEN

		WRITE(*,*) '!here!', SV

		ALLOCATE( A(rank, rank), B(N-rank,N-rank), C(rank,N-rank), &
			  Q(N-rank,N-rank), Ag(rank,rank), Qg(N-rank,N-rank) )

		A = M(1:rank,1:rank)
		B = M(rank+1:N,rank+1:N)
		C = M(1:rank,rank+1:N)

		CALL inverse( rank, lambda, A, Ag )

		Q = B - MATMUL( conjg( transpose( C ) ), MATMUL( Ag, C ) ) 
		CALL inverse( N-rank, lambda, Q, Qg )

		Mg(rank+1:N,1:rank) = -MATMUL( Qg, &
					       MATMUL( CONJG( TRANSPOSE( C ) ), Ag ) )
		Mg(1:rank,rank+1:N) = -MATMUL( Ag, &
					       MATMUL( C, Qg ) )
		Mg(1:rank,1:rank) = Ag + MATMUL( Ag, &
						 MATMUL( C, -Mg(rank+1:N,1:rank) ) )
		Mg(rank+1:N,rank+1:N) = Qg

	ELSE

		CALL inverse( N, lambda, M, Mg )

	END IF

END SUBROUTINE
