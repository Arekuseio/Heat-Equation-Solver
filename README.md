# Heat-Equation-Solver
The solution for one-dimensional heat equation using finite difference method (difference scheme with sigma weight). The approximation error is O(thau^2 + h^2), where thau is the time step and h is space step.

To solve the equation in the grid nodes, it is necessary to set the right side of the equation by the function f(x, t), as well as the boundary conditions of the problem for u(0, x) = mu, u(t, 0) = mu1, u(1, t) = mu2. The result is an approximate solution at the grid nodes and the solution time for the various grid steps.

At a fixed time step, a tridiagonal matrix is compiled to solve the equation at each internal node in space, the solution of such a system will be an approximate vector for solving the problem at a fixed time.

If you know the real solution, you can also calculate the maximum error among all the grid nodes.
