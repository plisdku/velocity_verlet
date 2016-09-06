function [xv, accelerations, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc)
% velocityVerlet3D     Calculate trajectory with driving force
%
% [xv, accel, ix_x, ix_y, ix_z, iv_x, iv_y, iv_z] = velocityVerlet3D(ts, xv0, accelFunc)
%
% accelFunc(t, x) returns the particle acceleration at one time and
% position, as a 3-vector [ax; ay; az].
%
% xv0 = [x0; y0; z0; vx0; vy0; vz0] are the initial conditions.
% ts should be uniformly-spaced but this is not checked.
%
% xv is the position and velocity at times ts(2:end).
% accel is the acceleration applied at times ts(1:end).
% Note that these are not the same length.
%
% accel should not be needed in practice but may be useful for debugging.
% (PCH uses it to compare values from accelFunc with trilinear
% interpolation done "by hand".)

Nt = length(ts) - 1;
dt = ts(2) - ts(1);

% System state: position and velocity
xs = zeros(Nt+1, 3);
vs = zeros(Nt+1, 3);
xs(1,1:3) = xv0(1:3);
vs(1,1:3) = xv0(4:6);

accelerations = zeros(Nt+1, 3);

currentAcceleration = accelFunc(ts(1), xs(1,1:3));
assert(isequal(size(currentAcceleration), [3, 1]));

accelerations(1,1:3) = currentAcceleration;

% Integrate!
for nn = 1:Nt
    
    xs(nn+1,:) = xs(nn,:) + vs(nn,:)*dt + 0.5*dt*dt*currentAcceleration';
    nextAcceleration = accelFunc(ts(nn+1), xs(nn+1,:));
    vs(nn+1,:) = vs(nn,:) + 0.5*dt*(currentAcceleration + nextAcceleration)';
    accelerations(nn+1,1:3) = nextAcceleration;
    
    currentAcceleration = nextAcceleration;
end

xv = [xs(:,1); xs(:,2); xs(:,3);
    vs(:,1); vs(:,2); vs(:,3)];

nn = 1:Nt+1;
ix_x = nn;
ix_y = ix_x + Nt + 1;
ix_z = ix_y + Nt + 1;
iv_x = ix_z + Nt + 1;
iv_y = iv_x + Nt + 1;
iv_z = iv_y + Nt + 1;



% i_x1 = @(nn) nn;
% i_x2 = @(nn) nn + Nt;
% i_x3 = @(nn) nn + 2*Nt;
% i_v1 = @(nn) nn + 3*Nt;
% i_v2 = @(nn) nn + 4*Nt;
% i_v3 = @(nn) nn + 5*Nt;

