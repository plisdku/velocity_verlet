function [G, dGdxv] = finalDistanceObjective3D(xDesired, xv)

Nt = length(xv)/6;

xFinal = xv([Nt, 2*Nt, 3*Nt]);  % TODO figure out best/least cumbersome indexing approach

displacement = xFinal - xDesired;
G = displacement'*displacement;
dGdxv = 0*xv;
dGdxv([Nt, 2*Nt, 3*Nt]) = 2*displacement;

