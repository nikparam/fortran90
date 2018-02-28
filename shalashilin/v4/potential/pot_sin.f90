SUBROUTINE potential_energy(x,params,V)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: V

	V = DSIN(x)
	RETURN

END SUBROUTINE potential_energy

SUBROUTINE diff_potential_energy(x,params,dV)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: dV

	dV = DCOS(x)
	RETURN

END SUBROUTINE diff_potential_energy

SUBROUTINE diff2_potential_energy(x,params,d2V)

	IMPLICIT NONE
	DOUBLE PRECISION, INTENT(IN) :: x
	DOUBLE PRECISION, INTENT(IN) :: params(15)
	DOUBLE PRECISION, INTENT(OUT) :: d2V

	d2V = -DSIN(x)
	RETURN

END SUBROUTINE diff2_potential_energy
