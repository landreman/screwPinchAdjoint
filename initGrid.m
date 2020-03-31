function grid = initGrid(N,params)

% This subroutine sets up a Chebyshev grid on [0,1] and rescales it to a
% grid for r.

% rho = r / a.

xgrid = ChebyshevGrid(N, 0, 1); % Arguments: # of grid points, xmin, xmax

r = xgrid.x * params.a;
ddr = params.a * xgrid.ddx;
d2dr2 = params.a * params.a * xgrid.d2dx2;

grid = struct('N',N,'rho',xgrid.x,'r',r,'ddrho',xgrid.ddx,'ddr',ddr,'d2drho2',xgrid.d2dx2,'d2dr2',d2dr2);

end