function [G, dGdxv] = sumObjective3D(xv)

G = sum(xv);
dGdxv = ones(size(xv));