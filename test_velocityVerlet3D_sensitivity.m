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
interpolant = @(E, x, y, z) interpn(x_grid, y_grid, z_grid, E, x, y, z, 'linear', 0);
accelFunc = @(Ex,Ey,Ez) @(t, xyz) [...
    interpolant(Ex, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ey, xyz(1), xyz(2), xyz(3)); ...
    interpolant(Ez, xyz(1), xyz(2), xyz(3)); ...
    ];

Nt = 100;
ts = linspace(0, 1, Nt);

[ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = get_Index3D(Nt);
ia_x = ix_x;
ia_y = ix_y;
ia_z = ix_z;

%%

xv0 = [0; 0; 0; .8; .499; 0];

[xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));

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
nn = 1:Nt;
accelInterpMatrix = sparse([nn, nn, nn, nn, nn, nn, nn, nn], ...
    [i000; i001; i010; i011; i100; i101; i110; i111], ...
    [w000; w001; w010; w011; w100; w101; w110; w111], ...
    Nt, numel(E_x));

ax = accelInterpMatrix * E_x(:);
ay = accelInterpMatrix * E_y(:);
az = accelInterpMatrix * E_z(:);

% To be sure we're right, compare the interpolated thing here to the accel
% variable that popped out of velocityVerlet3D.

checkClose(ax, accel(:,1));
checkClose(ay, accel(:,2));
checkClose(az, accel(:,3));
disp('Acceleration interpolation test PASSED!');

%% Compare matrix result to solver result

[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);

xv_mat = systemMatrix \ (-initMatrix*xv0 -accelMatrix*[ax; ay; az]);

[Gmat, DGmat] = sumObjective3D(xv_mat);

% Test a few dual sensitivities: an acceleration and the initial conds.
xv_dual = systemMatrix' \ DGmat;

Dxv0 = [1; 0; 0; 0; 0; 0];
dGdx0 = xv_dual' * (-initMatrix * Dxv0);

delta = 1e-6;
xv2 = systemMatrix \ (-initMatrix*(xv0 + delta*Dxv0) - accelMatrix*[ax;ay;az]);
[Gmat2, DGmat2] = sumObjective3D(xv2);

checkClose = @(a, b) assert(norm(a-b) < 1e-6);
checkClose(dGdx0, (Gmat2 - Gmat)/delta);
disp('Sum objective sensitivity to x0 test PASSED!')

%% Sensitivity with respect to ax(t)!  No interpolation yet.

[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
[xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));
[G, DG] = sumObjective3D(xv);

xv_dual = systemMatrix' \ DG;
dGda = xv_dual' * (-accelMatrix);

Da = zeros(3*Nt,1);
Da(5) = 1;
dGda5 = dGda * Da;

delta = 1e-6;
acceleration = [accel(:,1); accel(:,2); accel(:,3)];
xv2 = systemMatrix \ (-initMatrix*xv0 -accelMatrix*(acceleration + delta*Da));

G2 = sumObjective3D(xv2);

dGda5_meas = (G2-G)/delta;
checkClose(dGda5, dGda5_meas);
disp('Acceleration sensitivity test PASSED!');

%% Now at last the sensitivities with respect to E-fields!

[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
[xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));
[G, DG] = sumObjective3D(xv);

xv_dual = systemMatrix' \ DG;
dGda = xv_dual' * (-accelMatrix);

dGdax = dGda(ia_x);
dGday = dGda(ia_y);
dGdEx = reshape(dGdax * accelInterpMatrix, [Nx Ny Nz]);
dGdEy = reshape(dGday * accelInterpMatrix, [Nx Ny Nz]);

figure(3); clf
imagesc(x_grid, y_grid, dGdEx(:,:,1)')
axis xy image vis3d
colorbar
ax = axis;
hold on
plot(xv(ix_x), xv(ix_y), 'b-');
%quiver(xv(ix_x), xv(ix_y), dGdax', dGday', 'r'); % hard to interpret.

%% Check that accelInterpMatrix gives da/dE.



%% Obtain all the Ex gradients by exhaustion

ax1 = accel(:,1);
ay1 = accel(:,2);
az1 = accel(:,3);
%%

[systemMatrix, initMatrix, accelMatrix] = velocityVerletMatrices3D(ts);
[xv, accel] = velocityVerlet3D(ts, xv0, accelFunc(E_x, E_y, E_z));
acceleration = [accel(:,1); accel(:,2); accel(:,3)];

xv_mat = systemMatrix \ (-initMatrix*xv0 -accelMatrix*acceleration);
%%
delta = 1e-1;
dGdEx_meas = 0*E_x(:,:,1);

%for xx = 1:Nx
%for yy = 1:Ny
for xx = 5:8
for yy = 5:8
    fprintf('%i, %i\n', xx, yy);
    
    Ex2 = E_x;
    Ex2(xx,yy,1) = Ex2(xx,yy,1) + delta;
    [xv2, accel2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
    [G2] = sumObjective3D(xv2);
%     Ex2 = E_x;
%     Ex2(xx,yy,1) = Ex2(xx,yy,1) - delta;
%     [xv2] = velocityVerlet3D(ts, xv0, accelFunc(Ex2, E_y, E_z));
%     [G3] = sumObjective3D(xv2);
    
    dGdEx_meas(xx,yy) = (G2-G)/delta;
    
    ax2 = accelInterpMatrix * Ex2(:);
    ay2 = accelInterpMatrix * E_y(:);
    az2 = accelInterpMatrix * E_z(:);
    
    figure(1); clf
    plot(ax1);
    hold on
    plot(ax2);
    legend('Orig', 'Perturbed')
    pause

%     figure(10); clf
%     imagesc(x_grid, y_grid, Ex2(:,:,1)');
%     axis xy image vis3d
%     colorbar
%     ax = axis;
%     hold on
%     plot(xv2(ix_x), xv2(ix_y), 'b-');
%     pause(0.01)
    
    
end
end
%%
figure(1); clf
subplot(121)
imagesc(x_grid, y_grid, dGdEx(:,:,1)');
axis xy image
colorbar
title('Adjoint')
subplot(122)
imagesc(x_grid, y_grid, dGdEx_meas');
axis xy image
colorbar
title('Meas')

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











