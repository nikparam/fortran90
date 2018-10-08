program test

	implicit none

	integer :: N, i, j
	double precision :: lambda
	double complex, allocatable :: M(:,:), Mg(:,:), comp(:,:)

	N = 4
	lambda = 0.0
	allocate( M( N, N ), Mg( N, N ), comp(N,N) )
	DO i = 1, N
		M(i,i) = DCMPLX( TAN( 0.5 * ( i + j ) ), 0.0D0 )
		DO j = i+1, N
			M(i,j) = DCMPLX( COS( 0.5 * i ),  SIN( 0.25 * j ) )
			M(j,i) = conjg( M(i,j) )
		END DO
	END DO

	CALL block_inverse( N, lambda, M, Mg )

	comp = MATMUL( M, MATMUL( Mg, M ) )

	DO i = 1, N
		WRITE(*,*) M(i,:) - comp(i,:)
	END DO
	

end program test
