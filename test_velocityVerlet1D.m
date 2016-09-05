%% Paul's velocity Verlet with dual
%
% Goal: an EXACT adjoint gradient!
% Everything will be done in discrete fashion beginning with the velocity
% Verlet implementation.

checkClose = @(a, b) assert(norm(a-b) < 1e-6);

Nt = 10;
ts = linspace(0, 1, Nt+1);

[systemMatrix, initMatrix, accelMatrix, i_x, i_v, i_a] = velocityVerletMatrices(ts);

dt = ts(2)-ts(1);

%% Test case: initial velocity

acceleration = zeros(Nt+1, 1);
xv0 = [0; 1];

xv = systemMatrix \ (-initMatrix*xv0 - accelMatrix*acceleration);
xs = xv(i_x(1:Nt));
vs = xv(i_v(1:Nt));

figure(1); clf
plot(ts(2:end), xs, 'o-');
hold on
plot(ts(2:end), vs, 'rx');

checkClose(vs(end), 1.0);
checkClose(xs(end), 1.0);
disp('Initial velocity test PASSED');

%% Test case: constant acceleration

acceleration = ones(Nt+1, 1);
xv0 = [0; 0];

xv = systemMatrix \ (-initMatrix*xv0 - accelMatrix*acceleration);
xs = xv(i_x(1:Nt));
vs = xv(i_v(1:Nt));

figure(1); clf
plot(ts(2:end), xs, 'o-');
hold on
plot(ts(2:end), vs, 'rx');

checkClose(xs(end), 0.5);
checkClose(vs(end), 1.0);
disp('Constant acceleration test PASSED');

%% Test case: initial velocity, using iterative solver

accelFunc = @(t, x) 0.0;
xv0 = [0; 1];

[xv, ~, ~] = velocityVerlet1D(ts, xv0, accelFunc);
xs = xv(1:Nt);
vs = xv(Nt+1:end);

figure(1); clf
plot(ts(2:end), xs, 'o-');
hold on
plot(ts(2:end), vs, 'rx');

checkClose(xs(end), 1.0);
checkClose(vs(end), 1.0);
disp('Initial velocity with solver test PASSED');

%% Test case: constant acceleration, using iterative solver

accelFunc = @(t, x) 1.0;
xv0 = [0; 0];

[xv, ~, ~] = velocityVerlet1D(ts, xv0, accelFunc);
xs = xv(1:Nt);
vs = xv(Nt+1:end);

figure(1); clf
plot(ts(2:end), xs, 'o-');
hold on
plot(ts(2:end), vs, 'rx');

checkClose(xs(end), 0.5);
checkClose(vs(end), 1.0);
disp('Constant acceleration with solver test PASSED');

%% Test case: system sensitivity to initial conditions and applied accel.
% Use a simple quadratic objective function

accelFunc = @(t, x) 1.0;
xv0 = [0; 0];

[xv, ~, ~] = velocityVerlet1D(ts, xv0, accelFunc);
xs = xv(1:Nt);
vs = xv(Nt+1:end);

[G, DG] = finalDistanceObjective(1.0, xv);

% Test a few dual sensitivities: an acceleration and the initial conds.
delta = 1e-9;

xv_dual = systemMatrix' \ DG;

% Sensitivity to x0:

Dxv0 = [1; 0];
dGdx0 = xv_dual' * (-initMatrix * Dxv0);

xv2 = systemMatrix \ (-initMatrix*(xv0 + delta*Dxv0) - accelMatrix*acceleration);
G2 = finalDistanceObjective(1.0, xv2);
dGdx0_meas = (G2-G)/delta;

checkClose(dGdx0, dGdx0_meas);

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






