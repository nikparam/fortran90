SUBROUTINE change_var(NEQ,q,p,ksi,eta,omega,phase)

	INTEGER, INTENT(IN) :: NEQ
	DOUBLE PRECISION, INTENT(IN) :: q(NEQ), p(NEQ), omega(NEQ), phase(NEQ)
	DOUBLE COMPLEX, INTENT(OUT) :: ksi(NEQ), eta(NEQ)

	DO i=1,NEQ
		ksi(i) = DCMPLX( omega(i) * q(i), p(i) )
		eta(i) = DCMPLX( phase(i) - 0.5D0 * omega(i) * q(i) * q(i), -q(i) * p(i) )
!		eta(i) = DCMPLX( phase(i) - 0.5D0 * omega(i) * q(i) * q(i), 0.0 )
	END DO

	RETURN

END SUBROUTINE change_var
