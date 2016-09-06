%% Test trilinear weights

%% Single box, simple cases
% This is as simple as it gets.

xgrid = [0, 1];
ygrid = [0, 1];
zgrid = [0, 1];

% Test points: 000, 100, 010, 001, 111, center
x = [0, 1, 0, 0, .5];
y = [0, 0, 1, 0, .5];
z = [0, 0, 0, 1, .5];

[ix, iy, iz, w000, w001, w010, w011, w100, w101, w110, w111] = ...
    trilinear_weights(x, y, z, xgrid, ygrid, zgrid);

assert(all(ix == 1))
assert(all(iy == 1))
assert(all(iz == 1))

assert(all(w000 == [1 0 0 0 0.125]))
assert(all(w001 == [0 1 0 0 0.125]))
assert(all(w010 == [0 0 1 0 0.125]))
assert(all(w011 == [0 0 0 0 0.125]))
assert(all(w100 == [0 0 0 1 0.125]))
assert(all(w101 == [0 0 0 0 0.125]))
assert(all(w110 == [0 0 0 0 0.125]))
assert(all(w111 == [0 0 0 0 0.125]))
assert(all(w111 == [0 0 0 0 0.125]))

disp('Single box tests PASSED');

%% Eight boxes
% Pay attention.  :-)  I test weird edge cases in here.

xgrid = [0, 0.5, 1];
ygrid = [0, 0.5, 1];
zgrid = [0, 0.5, 1];

% Test points: 000, 100, 010, 001, 111, and center of 111 corner box
x = [0, 1, 0, 0, .75];
y = [0, 0, 1, 0, .75];
z = [0, 0, 0, 1, .75];

[ix, iy, iz, w000, w001, w010, w011, w100, w101, w110, w111] = ...
    trilinear_weights(x, y, z, xgrid, ygrid, zgrid);

% The fact that these tests pass is due to exactly accurate floating point
% division on nice fractions like 0.5... in general edge cases are risky.
assert(all(ix == [1, 2, 1, 1, 2]));
assert(all(iy == [1, 1, 2, 1, 2]));
assert(all(iz == [1, 1, 1, 2, 2]));

% To get these values, figure out first which box we're in and then what
% the corner weights ought to be.  This is actually a pretty strange test.
assert(all(w000 == [1 0 0 0 0.125]))
assert(all(w001 == [0 1 0 0 0.125]))
assert(all(w010 == [0 0 1 0 0.125]))
assert(all(w011 == [0 0 0 0 0.125]))
assert(all(w100 == [0 0 0 1 0.125]))
assert(all(w101 == [0 0 0 0 0.125]))
assert(all(w110 == [0 0 0 0 0.125]))
assert(all(w111 == [0 0 0 0 0.125]))
assert(all(w111 == [0 0 0 0 0.125]))

disp('Eight-box tests PASSED');
