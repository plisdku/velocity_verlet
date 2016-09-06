function [xv, i_x, i_v] = velocityVerlet1D(ts, xv0, accelFunc)
% velocityVerlet1D     Calculate trajectory with driving force
%
% [xv, i_x, i_v] = velocityVerlet1D(ts, xv0, accelFunc)
%
% accelFunc(t, x) returns the particle acceleration at one time and
% position.
%
% xv0 = [x0; v0] are the initial conditions.
% ts should be uniformly-spaced but this is not checked.

Nt = length(ts) - 1;
dt = ts(2) - ts(1);

% System state: position and velocity
xs = zeros(Nt+1, 1);
vs = zeros(Nt+1, 1);
xs(1) = xv0(1);
vs(1) = xv0(2);

currentAcceleration = accelFunc(ts(1), xs(1));

% Integrate!
for nn = 1:Nt
    
    xs(nn+1) = xs(nn) + vs(nn)*dt + 0.5*dt*dt*currentAcceleration;
    nextAcceleration = accelFunc(ts(nn+1), xs(nn+1));
    vs(nn+1) = vs(nn) + 0.5*dt*(currentAcceleration + nextAcceleration);
    
    currentAcceleration = nextAcceleration;
    
end

xv = [xs; vs];
i_x = 1:Nt+1;
i_v = i_x + Nt + 1;

