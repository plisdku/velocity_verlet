%% 

Nx          = 11;
Ny          = 11;
Nz          = 2;
x_grid = linspace(-1, 1, Nx);
y_grid = linspace(-1, 1, Ny);
z_grid = linspace(-1, 1, Nz);

%V           = ones(Nx, Ny, Nz);
%V           = 3*randn(Nx, Ny, Nz);
V           = repmat(peaks(Nx), [1 1 2]); assert(Nx == Ny);
V(:,:,2) = V(:,:,1);
E_x         = -centeredDiff(V, 1);
E_y         = -centeredDiff(V, 2);
E_z         = -centeredDiff(V, 3);

% put these out here so I can use them...
interp_E_x = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z, 'linear', 0);
interp_E_y = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z, 'linear', 0);
interp_E_z = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_z, x, y, z, 'linear', 0);

accelFunc = @(t, xyz) [...
    interp_E_x(xyz(1), xyz(2), xyz(3));
	interp_E_y(xyz(1), xyz(2), xyz(3));
	interp_E_z(xyz(1), xyz(2), xyz(3));
	];

Nt = 100;
ts = linspace(0, 1, Nt+1);

%%

xv0 = [0; 0; 0; .8; .499; 0];

[xv, accel, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc);

figure(1); clf
imagesc(x_grid, y_grid, V(:,:,1)');
axis xy image vis3d
ax = axis;
hold on
plot3(xv(ix_x), xv(ix_y), xv(ix_z), 'w.-');
quiver3(xv(ix_x), xv(ix_y), xv(ix_z), ...
    accel(:,1), accel(:,2), accel(:,3), 2, 'r')
axis(ax)

%% Now we get to figure out the trilinear interpolation part.
% Build a matrix to turn E-fields into accelerations along the path.

[iix, iiy, iiz, w000, w001, w010, w011, w100, w101, w110, w111] = ...
    trilinear_weights(xv(ix_x), xv(ix_y), xv(ix_z), ...
    x_grid, y_grid, z_grid);

% Turn this into a matrix to eat the E_x, E_y and E_z grids and produce
% the accelerations.
sz = [Nx Ny Nz];
i000 = sub2ind(sz, iix, iiy, iiz);
i001 = sub2ind(sz, iix+1, iiy, iiz);
i010 = sub2ind(sz, iix, iiy+1, iiz);
i011 = sub2ind(sz, iix+1, iiy+1, iiz);
i100 = sub2ind(sz, iix, iiy, iiz+1);
i101 = sub2ind(sz, iix+1, iiy, iiz+1);
i110 = sub2ind(sz, iix, iiy+1, iiz+1);
i111 = sub2ind(sz, iix+1, iiy+1, iiz+1);

% Test the accelerations: YES!  I am awesome.
%ax = E_x(i000).*w000 + E_x(i001).*w001 + E_x(i010).*w010 + E_x(i011).*w011 ...
%    + E_x(i100).*w100 + E_x(i101).*w101 + E_x(i110).*w110 + E_x(i111).*w111;

% This matrix should interpolate my E-field for me!  Nothing more needed.
nn = 1:Nt+1;
accelMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
    [i000; i001; i010; i011; i100; i101; i110; i111], ...
    [w000; w001; w010; w011; w100; w101; w110; w111], ...
    Nt+1, numel(E_x));

ax = accelMatrix * E_x(:);
ay = accelMatrix * E_y(:);
az = accelMatrix * E_z(:);

% To be sure we're right, compare the interpolated thing here to the accel
% variable that popped out of velocityVerlet3D.

%% Compare matrix result to solver result

[systemMatrix, initMatrix, accelMatrix, ...
    ix_x, ix_y, ix_z, ...
    iv_x, iv_y, iv_z, ...
    ia_x, ia_y, ia_z] = velocityVerletMatrices3D(ts);

xv_mat = systemMatrix \ (-initMatrix*xv0 -accelMatrix*[ax; ay; az]);

[Gmat, DGmat] = sumObjective3D(xv_mat);

% Test a few dual sensitivities: an acceleration and the initial conds.
delta = 1e-6;
xv_dual = systemMatrix' \ DGmat;

Dxv0 = [1; 0; 0; 0; 0; 0];
dGdx0 = xv_dual' * (-initMatrix * Dxv0);

xv2 = systemMatrix \ (-initMatrix*(xv0 + delta*Dxv0) - accelMatrix*[ax;ay;az]);
[Gmat2, DGmat2] = sumObjective3D(xv2);

disp('THIS IS WHAT WORKED')

checkClose = @(a, b) assert(norm(a-b) < 1e-6);
checkClose(dGdx0, (Gmat2 - Gmat)/delta);


%% Now at last the sensitivities with respect to E-fields!

[systemMatrix, initMatrix, accelMatrix, ...
    ix_x, ix_y, ix_z, ...
    iv_x, iv_y, iv_z, ...
    ia_x, ia_y, ia_z] = velocityVerletMatrices3D(ts);

[G, DG] = sumObjective3D(xv);

%nn = 1:Nt;
DG = DG([ix_x(nn), ix_y(nn), ix_z(nn), iv_x(nn), iv_y(nn), iv_z(nn)]);

%%

% Test a few dual sensitivities: an acceleration and the initial conds.
delta = 1e-3;

nn = 2:length(ts);
xv_dual = systemMatrix' \ DG;

%%
% Sensitivity to x0:

Dxv0 = [1; 0; 0; 0; 0; 0];
dGdx0 = xv_dual' * (-initMatrix * Dxv0);

xv2 = velocityVerlet3D(ts, xv0 + [delta; 0; 0; 0; 0; 0], accelFunc);
G2 = sumObjective3D(xv2);
dGdx0_meas = (G2-G)/delta;

checkClose(dGdx0, dGdx0_meas);
%%
% Sensitivity to v0:

Dxv0 = [0; 1];
dGdv0 = xv_dual' * (-initMatrix * Dxv0);

xv2 = systemMatrix \ (-initMatrix*(xv0 + delta*Dxv0) - accelMatrix*acceleration);
G2 = finalDistanceObjective(1.0, xv2);
dGdv0_meas = (G2-G)/delta;

checkClose(dGdx0, dGdv0_meas);

% Sensitivity to a(5):

Da = zeros(Nt+1, 1);
Da(5) = 1.0;
dGda = xv_dual' * (-accelMatrix * Da);

xv2 = systemMatrix \ (-initMatrix*xv0 - accelMatrix*(acceleration+delta*Da));
G2 = finalDistanceObjective(1.0, xv2);
dGda_meas = (G2-G)/delta;

checkClose(dGda, dGda_meas);

disp('Dual sensitivity tests PASSED');











