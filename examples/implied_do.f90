PROGRAM implied_do
	INTEGER :: i
	INTEGER, PARAMETER :: N = 5
!	READ(*,*) N
	
	REAL, DIMENSION(N) :: x
	x(1:N) = (/ (SQRT( REAL(i) ), i=1, N) /)
	DO i = 1, N
		WRITE(*,*) 'x(',i, ')= ', x(i)
	END DO
END PROGRAM implied_do
