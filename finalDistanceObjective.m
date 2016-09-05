function [G, dGdxv] = finalDistanceObjective(xDesired, xv)

Nt = length(xv)/2;

xFinal = xv(Nt);  % TODO figure out best/least cumbersome indexing approach

G = (xFinal - xDesired)^2;
dGdxv = 0*xv;
dGdxv(Nt) = 2*(xFinal - xDesired);

