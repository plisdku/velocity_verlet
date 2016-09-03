%% Wrapper

% Size of grid.  Set Nz = 2 for quasi-2D behavior.
% NOTE: removing Setup_Grid3D because it's harder to change the grid that
% way.
Nx          = 11;
Ny          = 11;
Nz          = 2;
x_grid = linspace(-1, 1, Nx);
y_grid = linspace(-1, 1, Ny);
z_grid = linspace(-1, 1, Nz);

V           = ones(Nx, Ny, Nz);
E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);

% "Boundary" in the time domain.
D_T         = [0, 1];
t            = 0:0.01:1;

nParticle   = 1;
x1_init     = 0.5*ones(1,nParticle);% 0.5 0.5 0.5];
x2_init     = 0*x1_init;% 0.4 0.3 0.2 0.1];
x3_init     = 0*x1_init;% 0.2 -0.4 -0.5 -0.6];
v1_init     = 0.2*ones(1,nParticle);% 0 0 0];
v2_init     = 0*x1_init;% 0.3 0.4 0.4 0.5];
v3_init     = v1_init;
x_v_init = Setup_Particle3D(x1_init, x2_init, x3_init, ...
    v1_init, v2_init, v3_init, nParticle);
m           = 1;
q           = 1;

x1_p         = -1*ones(1,nParticle);
x2_p         = 0*ones(1,nParticle);
x3_p         = 0*ones(1,nParticle);
objectiveFunction = @(x_v_fw) objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p);

%%

[dF_dE_x_dt1, dF_dE_y_d, dF_dE_z_d, G, D_G, x_v_fw] = dFdE_VV( ...
    x_grid, y_grid, z_grid, ...
    objectiveFunction, ...
    nParticle, t, ...
    x_v_init, ...
    E_x, E_y, E_z, ...
    m, q);

%%
[ix,iy,iz,~,~,~] = get_Index3D(nParticle);

figure(1); clf
plot(x_v_fw(:,ix), x_v_fw(:,iy), 'o-')

%% Make sure the function runs
% ... it doesn't.

% dFdV_VV( ...
%     x_grid, y_grid, z_grid, ...
%     x1_p, x2_p, x3_p, ...
%     nParticle, t, ...
%     x1_init, x2_init, x3_init, v1_init, v2_init, v3_init, ...
%     V, m, q);

%%
