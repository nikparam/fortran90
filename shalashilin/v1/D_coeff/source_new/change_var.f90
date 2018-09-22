SUBROUTINE change_var(NumG,q,p,ksi,eta,omega,phase)

	IMPLICIT NONE

	INTEGER, INTENT(IN) :: NumG
	DOUBLE PRECISION, INTENT(IN) :: q(NumG), p(NumG), omega(NumG), phase(NumG)
	DOUBLE COMPLEX, INTENT(OUT) :: ksi(NumG), eta(NumG)

	INTEGER :: i

	DO i=1,NumG
		ksi(i) = DCMPLX( omega(i) * q(i), p(i) )
		eta(i) = DCMPLX( phase(i) - 0.5D0 * omega(i) * q(i) * q(i), 0.0 )
	END DO

	RETURN

END SUBROUTINE change_var
