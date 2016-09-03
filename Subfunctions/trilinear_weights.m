function [i_x, i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111] = trilinear_weights(x,y,z, x_grid, y_grid, z_grid, d_x, d_y, d_z, nParticle)
% trilinear_weights    Indices and weights for trilinear interpolation
%
% trilinear_weights(x, y, z, xgrid, ygrid, zgrid, dx, dy, dz, nParticle)
% Calculate weights and indices needed to interpolate a function defined
% on a regular 3D grid to nParticle lists of coordinates.
%
% Inputs:
%    x, y, z            N-element 
% 
% Outputs:
%   i_x, i_y, i_z       Indices to evaluate function at
%   w000, ..., w111     Weights to assign to function values
%
% TODO: Remove multi-particle functionality (separation of concerns)
% TODO: Infer dx, dy, dz from x_grid, y_grid, z_grid (simplification)

% TODO: what does "log" mean in "log_vec"?
log_vec = ((x >= x_grid(1)) & (x <= x_grid(end)) & (y >= y_grid(1)) & (y <= y_grid(end)) & (z >= z_grid(1)) & (z <= z_grid(end)));

i_x = ones(size(x,2),nParticle);
i_y = ones(size(y,2),nParticle);
i_z = ones(size(z,2),nParticle);

if nParticle == 1
    
    i_x(log_vec) = floor( (x(log_vec) - x_grid(1))/d_x) +1;
    i_y(log_vec) = floor( (y(log_vec) - y_grid(1))/d_y) +1;
    i_z(log_vec) = floor( (z(log_vec) - z_grid(1))/d_z) +1;

    x_c = x_grid(i_x);
    y_c = y_grid(i_y);
    z_c = z_grid(i_z);
    
else

    for iParticle = 1:nParticle
        i_x(log_vec(iParticle,:),iParticle) = floor( (x(iParticle,log_vec(iParticle,:)) - x_grid(1)) / d_x ) + 1;
        i_y(log_vec(iParticle,:),iParticle) = floor( (y(iParticle,log_vec(iParticle,:)) - y_grid(1)) / d_y ) + 1;
        i_z(log_vec(iParticle,:),iParticle) = floor( (z(iParticle,log_vec(iParticle,:)) - z_grid(1)) / d_z ) + 1; 
    end
    % get the Value of the corner

    x_c = x_grid(i_x);
    y_c = y_grid(i_y);
    z_c = z_grid(i_z);

end

% Recover weights

wx = ( x - x_c' ) ./ d_x;
wy = ( y - y_c' ) ./ d_y;
wz = ( z - z_c' ) ./ d_z; 

w000 = (1-wx).*(1-wy).*(1-wz);
w100 = wx.*(1-wy).*(1-wz);
w010 = (1-wx).*wy.*(1-wz);
w110 = wx.*wy.*(1-wz);
w001 = (1-wx).*(1-wy).*wz;
w101 = wx.*(1-wy).*wz;
w011 = (1-wx).*wy.*wz;
w111 = wx.*wy.*wz;

w000(~log_vec) = 0;
w100(~log_vec) = 0;
w010(~log_vec) = 0;
w110(~log_vec) = 0;

w001(~log_vec) = 0;
w101(~log_vec) = 0;
w011(~log_vec) = 0;
w111(~log_vec) = 0;

end