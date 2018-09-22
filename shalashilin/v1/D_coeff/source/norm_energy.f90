SUBROUTINE norm_energy( NEQ, q, p, D, N, E, LE, omega, params, phase, npts, x, wts, dE )

	DOUBLE PRECISION, INTENT(IN) :: q(NEQ), p(NEQ), omega(NEQ), params(15), phase(NEQ), dE
	DOUBLE COMPLEX, INTENT(IN) :: D(NEQ)
	DOUBLE PRECISION, INTENT(OUT) :: N, E, LE
	DOUBLE COMPLEX :: ksi(NEQ), eta(NEQ), S(NEQ,NEQ), H(NEQ,NEQ), M1(NEQ,NEQ), L(NEQ,NEQ) 

	CALL change_var( NEQ, q, p, ksi, eta, omega, phase )
	CALL overlap( NEQ, ksi, eta, omega, S )
	CALL hamiltonian( NEQ, ksi, eta, omega, params, npts, x, wts, dE, S, H, L, M1 )

	N = DOT_PRODUCT(D, MATMUL(S,D))
	E = DOT_PRODUCT(D, MATMUL(H,D))/N
	LE = DOT_PRODUCT(D, MATMUL(L,D))/N

	RETURN

END SUBROUTINE
