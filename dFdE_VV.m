%% Function VV

function [dF_dE_x_d, dF_dE_y_d, dF_dE_z_d, G, D_G, x_v_fw] = dFdE_VV(...
    x_grid, y_grid, z_grid, ...
    objectiveFunction, nParticle, t, x_v_init, E_x, E_y, E_z, m, q)
%% Init
%[x_grid, y_grid, z_grid, d_x, d_y, d_z]     = Setup_Grid3D(N);

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]  = get_Index3D(nParticle);

interp_E_x = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z,...
    'linear', 0);
interp_E_y = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z,...
    'linear', 0);
interp_E_z = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_z, x, y, z,...
    'linear', 0);

delta_t = t(2) - t(1);

%% Solve Systems
[x_v_fw, S] = VV3D(q, m, t, delta_t, nParticle, interp_E_x, ...
    interp_E_y, interp_E_z, x_v_init);

[G, D_G] = objectiveFunction(x_v_fw);
%[G, D_G] = objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p);

[xd] = dualVV(nParticle, D_G, S, t);

%% Interpolation in Time and Weight Recovery
%
x_d = spline(t, x_v_fw(:,id_x1)', t);
y_d = spline(t, x_v_fw(:,id_x2)', t);
z_d = spline(t, x_v_fw(:,id_x3)', t);

[i_x, i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111] = ...
    trilinear_weights(x_d' ,y_d' ,z_d', x_grid, y_grid, z_grid, nParticle);

%% Get the Value of Integrant
%
d_t_s       = -centeredDiff(t)';
dF_dE_x_t   = bsxfun(@times,d_t_s, xd(:,id_v1)*q/m);
dF_dE_y_t   = bsxfun(@times,d_t_s, xd(:,id_v2)*q/m);
dF_dE_z_t   = bsxfun(@times,d_t_s, xd(:,id_v3)*q/m);

%% Calculate Sensitivity with Dual Function 
%
[dF_dE_x_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_x_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );
[dF_dE_y_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_y_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );
[dF_dE_z_d]  = multipleRestriction3D(i_x, i_y, i_z,...
    w000, w100, w010,w110, w001, w101, w011, w111, ...
    dF_dE_z_t, size(E_x,1), size(E_x,2), size(E_x,3), nParticle );

end