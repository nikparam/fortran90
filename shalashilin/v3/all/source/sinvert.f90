SUBROUTINE inverse(NumG, lambda, S, SI)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: lambda
	DOUBLE COMPLEX, INTENT(IN) :: S(NumG,NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: SI(NumG,NumG)

	INTEGER :: i, LW, LIW, LRW, astatus, INFO
	DOUBLE PRECISION :: Sd(NumG)
	DOUBLE COMPLEX :: TEMP(NumG,NumG), Sdm(NumG,NumG), U(NumG,NumG), VT(NumG,NumG), &
			  TMP(NumG,NumG)

	INTEGER, ALLOCATABLE :: IWORK(:)
	DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
	DOUBLE COMPLEX, ALLOCATABLE :: WORK(:)


	LW = NumG*( NumG + 3 )
	LIW = 8 * NumG
	LRW = NumG * ( 5 * NumG + 7 )
	ALLOCATE(WORK(LW), IWORK(LIW), RWORK(LRW), STAT=astatus)

	TEMP(:,:) = S(:,:)
	CALL ZGESDD( 'A', NumG, NumG, TEMP, NumG, Sd, U, NumG, VT, NumG, WORK, LW, RWORK, IWORK, INFO )

	IF ( INFO .NE. 0 ) THEN
		WRITE(*,*) 'Something went wrong'
		WRITE(*,*) '____________________'
		DO i = 1, NumG
			WRITE(*,*) S(i,:)
		END DO
		WRITE(*,*) '____________________'
	END IF 

	Sdm(1:NumG,1:NumG) = (0.0D0, 0.0D0)
	DO i = 1, NumG
		Sdm(i,i) = Sd(i) / ( Sd(i) * Sd(i) + lambda**2 )
!		Sdm(i,i) = 1.0D0 / Sd(i)
	END DO

	SI = MATMUL( TRANSPOSE( CONJG( VT ) ), MATMUL( Sdm, TRANSPOSE( CONJG( U ) ) ) )

	RETURN
END SUBROUTINE inverse
