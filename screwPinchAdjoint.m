params = struct('a',1,'R0',1000,'psi0',5/(2*pi),'mu0',4*pi*(1e-7));

Ns = [17, 33, 65];

figure(11)
clf
numRows = 2;
numCols = 2;

for j = 1:numel(Ns)
    N = Ns(j);
    grid = initGrid(N, params);

%{
alpha = 1.4;

f = sin(alpha*grid.x);
dfdx_analytic = alpha*cos(alpha*grid.x)

dfdx_diffMat = grid.ddx * f

difference = dfdx_analytic - dfdx_diffMat

d2fdx2_analytic = -f*alpha*alpha

d2fdx2_diffMat = grid.d2dx2 * f

difference = d2fdx2_analytic - d2fdx2_diffMat
%}

    psi = solveEquilibrium(params, grid);
    
    [xi_psi_forward, xi_psi_adjoint] = solve_m0_n0(params, grid, psi);
    
    subplot(numRows,numCols,1)
    plot(grid.r, psi,'.-','displayname',['N=',num2str(N)])
    hold on

    subplot(numRows,numCols,2)
    plot(psi, xi_psi_forward,'.-','displayname',['N=',num2str(N)])
    hold on

    subplot(numRows,numCols,3)
    plot(psi, xi_psi_adjoint,'.-','displayname',['N=',num2str(N)])
    hold on
end

subplot(numRows,numCols,1)
xlabel('r')
ylabel('\psi')
title('Unperturbed Equilibrium')
legend show
set(legend,'location','southeast')

subplot(numRows,numCols,2)
xlabel('\psi')
ylabel('\xi^\psi')
title('Reproducing Figure 6.1')
legend show
set(legend,'location','southeast')

subplot(numRows,numCols,3)
xlabel('\psi')
ylabel('\xi^\psi')
title('Reproducing Figure 6.2')
legend show
set(legend,'location','southeast')
%xlim([0,0.8])
%ylim([0,
