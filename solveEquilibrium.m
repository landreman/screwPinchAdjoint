function psi = solveEquilibrium(params, grid)

% Initial guess is a quadratic function:
psi = params.psi0 * grid.rho .* grid.rho;

options = optimoptions('fsolve','Display','iter','FiniteDifferenceType','central','StepTolerance',1e-16,'FunctionTolerance',1e-16);
psi = fsolve(@residual, psi, options)

    function res = residual(psi)
        % Handle the interior:
        psiPrime = grid.ddr * psi;
        iota = iota_profile(psi / params.psi0);
        res = grid.ddr * (params.mu0 * p_profile(psi / params.psi0) ) ...
            + 1 ./ (2 * grid.r .* grid.r) .* (grid.ddr * (psiPrime .^ 2)) ...
            - (psiPrime .^ 2) ./ (grid.r .* grid.r .* grid.r) ...
            + iota .* psiPrime ./ (grid.r * params.R0 * params.R0) .* (grid.ddr * (grid.r .* iota .* psiPrime));
        
        % Handle the boundary conditions:
        res(1) = psi(1);
        res(end) = psi(end) - params.psi0;
        %res
    end

fprintf('Residuals at final psi:\n')
residual(psi)

end