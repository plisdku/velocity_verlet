%% Test 3D velocity verlet solver

checkClose = @(a, b) assert(norm(a-b) < 1e-6);

Nt = 10;
ts = linspace(0, 1, Nt+1);

[systemMatrix, initMatrix, accelMatrix,...
    ix_x, ix_y, ix_z, ...
    iv_x, iv_y, iv_z, ...
    ia_x, ia_y, ia_z] = velocityVerletMatrices3D(ts);

dt = ts(2)-ts(1);

%% Test case: initial velocity

acceleration = zeros(3*(Nt+1), 1);
xv0 = [0; 0; 0; 1; 2; 3];

xv = systemMatrix \ (-initMatrix*xv0 - accelMatrix*acceleration);

% figure(1); clf
% plot(ts(2:end), xv(ix_x), 'o-');
% hold on
% plot(ts(2:end), xv(iv_x), 'rx');

checkClose(xv(iv_x(end)), 1.0);
checkClose(xv(iv_y(end)), 2.0);
checkClose(xv(iv_z(end)), 3.0);
checkClose(xv(ix_x(end)), 1.0);
checkClose(xv(ix_y(end)), 2.0);
checkClose(xv(ix_z(end)), 3.0);
disp('Initial velocity test PASSED');

%% Test case: constant acceleration

unos = ones(Nt+1,1);
acceleration = [unos; 2*unos; 3*unos];

xv0 = [0; 0; 0; 0; 0; 0];

xv = systemMatrix \ (-initMatrix*xv0 - accelMatrix*acceleration);

% figure(1); clf
% plot(ts(2:end), xv(ix_x), 'o-');
% hold on
% plot(ts(2:end), xv(iv_x), 'rx');

checkClose(xv(ix_x(end)), 0.5);
checkClose(xv(iv_x(end)), 1.0);
checkClose(xv(ix_y(end)), 1.0);
checkClose(xv(iv_y(end)), 2.0);
checkClose(xv(ix_z(end)), 1.5);
checkClose(xv(iv_z(end)), 3.0);
disp('Constant acceleration test PASSED');

%% Test case: initial velocity, using iterative solver

accelFunc = @(t, x) [0; 0; 0];
xv0 = [0; 0; 0; 1; 2; 3];

[xv, ~, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc);

% figure(1); clf
% plot(ts(2:end), xv(i_x1(1:Nt)), 'o-');
% hold on
% plot(ts(2:end), xv(i_v1(1:Nt)), 'rx');

checkClose(xv(ix_x(end)), 1.0);
checkClose(xv(iv_x(end)), 1.0);
checkClose(xv(ix_y(end)), 2.0);
checkClose(xv(iv_y(end)), 2.0);
checkClose(xv(ix_z(end)), 3.0);
checkClose(xv(iv_z(end)), 3.0);
disp('Initial velocity with solver test PASSED');

%% Test case: constant acceleration, using iterative solver

accelFunc = @(t, x) [1.0; 2.0; 3.0];
xv0 = [0; 0; 0; 0; 0; 0];

[xv, ~, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc);

% figure(1); clf
% plot(ts(2:end), xv(i_x1(1:Nt)), 'o-');
% hold on
% plot(ts(2:end), xv(i_v1(1:Nt)), 'rx');

checkClose(xv(ix_x(end)), 0.5);
checkClose(xv(iv_x(end)), 1.0);
checkClose(xv(ix_y(end)), 1.0);
checkClose(xv(iv_y(end)), 2.0);
checkClose(xv(ix_z(end)), 1.5);
checkClose(xv(iv_z(end)), 3.0);
disp('Constant acceleration with solver test PASSED');




