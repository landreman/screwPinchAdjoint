function [xi_psi_forward, xi_psi_adjoint] = solve_m0_n0(params, grid, psi)

% Define some shorthand:
R0 = params.R0;
r = grid.r;

d_psi_d_r = grid.ddr * psi;

d2_psi_d_r2 = grid.d2dr2 * psi;

iota = iota_profile(psi / params.psi0);

p = p_profile(psi / params.psi0);

d_iota_d_r = grid.ddr * iota;

factor = R0 * R0 + r .* r .* iota .* iota;

delta_F = - grid.ddr * p;

B1 = (R0 * R0 - r .* r .* iota .* (iota + 2 * r .* d_iota_d_r)) ./ (r .* factor) ...
    - 2 * d2_psi_d_r2 ./ d_psi_d_r;

B2 = ((3 * R0 * R0 - r .* r .* iota .* iota) .* d_psi_d_r - 2 * r .* R0 .* R0 .* d2_psi_d_r2) ...
    ./ ( r .* r .* factor .* d_psi_d_r);

B3 = - params.mu0 * r .* r .* delta_F ./ ((1 + r .* r .* iota .* iota / (R0 * R0)) .* d_psi_d_r .* d_psi_d_r);

% Form the interior of the matrix:
matrix = grid.d2dr2 - diag(B1) * grid.ddr - diag(B2);

% Handle boundary conditions:
matrix(1,:) = 0;
matrix(1,1) = 1;
matrix(end,:) = 0;
matrix(end,end) = 1;

% Solve the forward problem:
rhs = zeros(grid.N,1);
rhs(end) = 1;

xi_r = matrix \ rhs;
xi_psi_forward = xi_r .* d_psi_d_r;

% Solve the adjoint problem:
rhs = B3;
rhs(1) = 0;
rhs(end) = 0;

xi_r = matrix \ rhs;
xi_psi_adjoint = xi_r .* d_psi_d_r;

end