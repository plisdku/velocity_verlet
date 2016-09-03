%% Function VV

function [dF_dE_x_fw, dF_dE_x_fw2, G, D_G, x_v_fw] = dFdE_VV_fw(x_grid, y_grid, z_grid, x1_p, x2_p, x3_p, nParticle, t, x1_init, x2_init, x3_init, v1_init, v2_init, v3_init, E_x, E_y, E_z, m, q)
%% Init
%[x_grid, y_grid, z_grid, d_x, d_y, d_z]     = Setup_Grid3D(N);

[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]  = get_Index3D(nParticle);

[x_v_init]                                  = Setup_Particle3D(...
    x1_init,x2_init, x3_init, v1_init, v2_init, v3_init, nParticle);

interp_E_x = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z,...
    'linear', 0);
interp_E_y = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z,...
    'linear', 0);
interp_E_z = @(x, y, z) interpn(x_grid, y_grid, z_grid, E_z, x, y, z,...
    'linear', 0);

delta_t = t(2) - t(1);

D_T = [t(1) t(end)];

%% Solve System

[x_v_fw, S] = VV3D(q, m, t, delta_t, nParticle, interp_E_x, ...
    interp_E_y, interp_E_z, x1_init, x2_init, x3_init, v1_init,...
    v2_init, v3_init);


[t_fw,x_v_fw_2] = forwardSystem_3D(interp_E_x, interp_E_y, interp_E_z,...
                x_v_init, D_T, nParticle, q, m);

[G, D_G] = objectiveFunction_3D(x_v_fw, nParticle, x1_p, x2_p, x3_p);

dF_dE_x_fw = zeros(length(x_grid),length(y_grid),length(z_grid));
dF_dE_x_fw2 = zeros(length(x_grid),length(y_grid),length(z_grid));
G_t = G;

G = sum(G);

u = 1;
for ii = 1:N%17:23%1:length(x_grid)
    ii
    for jj = 1:N%13:N%1:length(y_grid)
        %jj
        for kk = 1:N%14:18%1:length(z_grid)
            %kk
            eps = 1e-9;
            
            
            E_x_pert = E_x;
            
            E_x_pert(ii,jj,kk) =  E_x_pert(ii,jj,kk) - 0.5 * eps;
            
            interp_E_x_pert = @(x, y, z) interpn(x_grid, y_grid, z_grid,...
                E_x_pert, x, y, z, 'linear', 0);
            
%             [t_Fw2, x_v_fw_pert] = forwardSystem_3D(interp_E_x_pert, ...
%                 interp_E_y, interp_E_z, x_v_init, D_T, nParticle, q, m);
            
            [x_v_fw_pert, ~] = VV3D(q, m, t, delta_t, nParticle, interp_E_x_pert, ...
                interp_E_y, interp_E_z, x1_init, x2_init, x3_init, v1_init,...
                v2_init, v3_init);
            
            [t_fw,x_v_fw_pert2] = forwardSystem_3D(interp_E_x_pert, interp_E_y, interp_E_z,...
                x_v_init, D_T, nParticle, q, m);
            
            [G_pert, D_G_pert] = objectiveFunction_3D(x_v_fw_pert, nParticle, x1_p, x2_p, x3_p);
            %G_pert_test1(ii+(jj-1)+(kk-1)) = sum(G_pert(id_x1));
            G_pert3 = sum(G_pert);
            
            [G_pert2, D_G_pert2] = objectiveFunction_3D(x_v_fw_pert2, nParticle, x1_p, x2_p, x3_p);
            %G_pert_test2(ii+(jj-1)+(kk-1)) = sum(G_pert2(id_x1));
%             G_pert_test(u) = G_pert3 - sum(G_pert2(id_x1));
%             u = u+1;
            G_pert_2 = sum(G_pert2);
            
            G_pert - G_pert2
            
            
           
            
            
            E_x_pert2 = E_x;
            E_x_pert2(ii,jj,kk) = E_x(ii,jj,kk) + 0.5*eps;
            
            interp_E_x_pert2 = @(x, y, z) interpn(x_grid, y_grid, z_grid,...
                E_x_pert2, x, y, z, 'linear', 0);
            
            
            
            [x_v_fw_pert_m2, ~] = VV3D(q, m, t, delta_t, nParticle, interp_E_x_pert2, ...
                interp_E_y, interp_E_z, x1_init, x2_init, x3_init, v1_init,...
                v2_init, v3_init);
            
            [t_fw,x_v_fw_pert2_m2] = forwardSystem_3D(interp_E_x_pert2, interp_E_y, interp_E_z,...
                x_v_init, D_T, nParticle, q, m);
            
            
            
            
            [G_pert_m2, D_G_pert] = objectiveFunction_3D(x_v_fw_pert_m2, nParticle, x1_p, x2_p, x3_p);
            %G_pert_test1(ii+(jj-1)+(kk-1)) = sum(G_pert(id_x1));
            G_pert3_m2 = sum(G_pert_m2);
            
            [G_pert2_m2, D_G_pert2] = objectiveFunction_3D(x_v_fw_pert2_m2, nParticle, x1_p, x2_p, x3_p);
            
            
            
            
            dF_dE_x_fw2(ii,jj,kk) = (sum(G_pert2) - sum(G_pert2_m2)) / eps;

            
            dF_dE_x_fw(ii,jj,kk) = (G_pert3 - G_pert3_m2) / eps;

            
        end
        
    end
end





end