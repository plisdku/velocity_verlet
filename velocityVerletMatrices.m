function [systemMatrix, A0, B, i_x, i_v, i_a] = velocityVerletMatrices(ts)
% velocityVerlet1D     Matrices and indexing for 1D velocity Verlet method
%
% systemMatrix, initialMatrix, accelMatrix, i_x, i_v, i_a] = velocityVerlet1D(ts)
%
% Let xv be the position and velocity vector we seek.  Then
%
% systemMatrix*xv + initialMatrix*[x0; v0] + accelMatrix*accelerations
%
% can be solved for xv.  Here,
%   length(ts) = Nt + 1
%   size(xv) = [2*Nt, 1]
%   size(accelerations) = [Nt+1, 1]     i.e. [a0; a1; ...; aNt]
%
% Use the indexing functions to get out the positions and velocities, i.e.
%   x(n) = xv(i_x(n))
%   v(n) = xv(i_v(n))
%

Nt = length(ts)-1;
dt = ts(2) - ts(1);

%% A single timestep:
%
% [x(n+1); v(n+1)] = A*[x(n); v(n)] + B*[a(n); a(n+1)]
%
% A = [1, dt; 0, 1]
% B = [0.5*dt^2, 0; 0.5*dt, 0.5*dt]
%
% x is position, v is velocity, and a is acceleration.
% We can build the difference equation for all timesteps out of these
% pieces, and build the 3D version out of that.  Actually since the three
% directions are uncoupled in the primal system we probably don't need to
% ever build the full system matrix for x, y and z.

% Indexing functions
i_x = @(timestep) timestep;
i_v = @(timestep) timestep + Nt;
i_a = @(timestep) timestep;

%% System matrix for 1D system
% We probably won't need to turn it into a matrix for a 3D system, but we
% could repeat it in block diagonal fashion.
%
% The system matrix updates timestep 0 to timestep 1, and so on, finally
% reaching timestep Nt.  The initial conditions are part of the right-hand
% side along with the driving forces at timesteps 0 through Nt.  NOTE that
% there are Nt+1 accelerations.

% ------- A
% A = [1, dt; 0, 1]

unos = @(n) ones(1,n);

ii_A = [i_x(2:Nt), i_x(2:Nt), i_v(2:Nt)];
jj_A = [i_x(1:Nt-1), i_v(1:Nt-1), i_v(1:Nt-1)];
vv_A = [unos(Nt-1), dt*unos(Nt-1), unos(Nt-1)];

% ------- Diagonal
% It's a big old identity matrix.

ii_diag = [i_x(1:Nt), i_v(1:Nt)];
jj_diag = [i_x(1:Nt), i_v(1:Nt)];
vv_diag = [-unos(Nt), -unos(Nt)];

% ------- Build the system matrix!

numRows = 2*Nt;
numCols = 2*Nt;

systemMatrix = sparse([ii_A, ii_diag], [jj_A, jj_diag], [vv_A, vv_diag],...
    numRows, numCols);
%spy(systemMatrix)

%% Right-hand side

% ------- B
% B = [0.5*dt^2, 0; 0.5*dt, 0.5*dt]
% The B matrix is 2*Nt x Nt+1.  Yeah it's weird-shaped.  It turns Nt+1
% accelerations into Nt timestep updates.

ii_B = [i_x(1:Nt), i_v(1:Nt), i_v(1:Nt)];
jj_B = [i_a(1:Nt), i_a(1:Nt), i_a(2:Nt+1)];
vv_B = [0.5*dt*dt*unos(Nt), 0.5*dt*unos(Nt), 0.5*dt*unos(Nt)];

numRows = 2*Nt;
numCols = Nt+1;
B = sparse(ii_B, jj_B, vv_B, numRows, numCols);
%figure(3); clf
%imagesc(B)

% ------- Initial conditions
% The first timestep needs to take initial conditions:
% [x(1); v(1)] = A*[x(0); v(0)] + B*[a(0); a(1)]
% So we need an RHS vector of the right size, 2*Nt x 1.
% Let's get RHS = A0 * [x(0); v(0)], with A0 being 2*Nt x 2.

numRows = 2*Nt;
numCols = 2;

% A = [1, dt; 0, 1]
ii_A0 = [i_x(1), i_x(1), i_v(1)];
jj_A0 = [1, 2, 2];
vv_A0 = [1, dt, 1];

A0 = sparse(ii_A0, jj_A0, vv_A0, numRows, numCols);

