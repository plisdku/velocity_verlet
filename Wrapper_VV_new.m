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
V           = 3*randn(Nx, Ny, Nz);
V(:,:,2) = V(:,:,1);
E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);

% put these out here so I can use them...
interp_E_x = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z, 'linear', 0);
interp_E_y = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z, 'linear', 0);
interp_E_z = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_z, x, y, z, 'linear', 0);

% "Boundary" in the time domain.
D_T         = [0, 1];
t            = 0:0.01:1;

nParticle   = 1;
x1_init     = -0.9*ones(1,nParticle);% 0.5 0.5 0.5];
x2_init     = 0*x1_init;% 0.4 0.3 0.2 0.1];
x3_init     = 0*x1_init;% 0.2 -0.4 -0.5 -0.6];
v1_init     = 2*ones(1,nParticle);% 0 0 0];
v2_init     = 0*x1_init;% 0.3 0.4 0.4 0.5];
v3_init     = 0*v1_init;
x_v_init = Setup_Particle3D(x1_init, x2_init, x3_init, ...
    v1_init, v2_init, v3_init, nParticle);
m           = 1;
q           = 1;

x1_p         = -1*ones(1,nParticle);
x2_p         = 0*ones(1,nParticle);
x3_p         = 0*ones(1,nParticle);
objectiveFunction = @(x_v_fw) objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p);

%%

[dF_dE_x_d, dF_dE_y_d, dF_dE_z_d, G, D_G, x_v_fw] = dFdE_VV( ...
    x_grid, y_grid, z_grid, ...
    objectiveFunction, ...
    nParticle, t, ...
    x_v_init, ...
    E_x, E_y, E_z, ...
    m, q);

%% Plot the E-field sensitivities...

figure(1); clf
imagesc(x_grid, y_grid, dF_dE_x_d(:,:,1)');
axis xy image
colorbar
xlabel('x')
ylabel('y')

%%
[ix,iy,iz,~,~,~] = get_Index3D(nParticle);

xs = x_v_fw(:,ix);
ys = x_v_fw(:,iy);
zs = x_v_fw(:,iz);

myEx = interp_E_x(xs, ys, zs);
myEy = interp_E_y(xs, ys, zs);

figure(2); clf
plot(xs, ys, 'o-')
xlim([-1,1])
ylim([-1,1])
%axis xy image
hold on
quiver(xs, ys, myEx, myEy)
xlabel('x')
ylabel('y')
title('Trajectory through random voltage')
legend('Particle', 'E-field')

%figure(2); clf
%plot3(xs, ys, zs, 'o-')

%% Make sure the function runs
% ... it doesn't.

% dFdV_VV( ...
%     x_grid, y_grid, z_grid, ...
%     x1_p, x2_p, x3_p, ...
%     nParticle, t, ...
%     x1_init, x2_init, x3_init, v1_init, v2_init, v3_init, ...
%     V, m, q);

%%
