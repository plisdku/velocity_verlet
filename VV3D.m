function [x_fw, S, x_v_fwold] = VV3D(q, m, t, delta_t, nParticle, interp_E_x, interp_E_y, interp_E_z, x1_init, x2_init, x3_init, v1_init, v2_init, v3_init)

[x_v_init] = Setup_Particle3D(x1_init,x2_init, x3_init, v1_init, v2_init, v3_init, nParticle);
[id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]  = get_Index3D(nParticle);


% x_grid = linspace(-1,1,11);
% y_grid = x_grid;
% z_grid = x_grid;
%E_x   = ones(11,11,11)*2;
%E_y   = E_x;
%E_y   = zeros(11,11,11);
%E_z = E_y;

%interp_E_x = @(x,y,z) interpn(x_grid, y_grid, z_grid, E_x, x, y, z, 'linear', 0);
%interp_E_y = @(x,y,z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z, 'linear', 0);
%interp_E_z = @(x,y,z) interpn(x_grid, y_grid, z_grid, E_y, x, y, z, 'linear', 0);



%xv_fw = zeros(length(t)*2,1);

x1_fw = zeros(length(t),nParticle);
v1_fw = zeros(length(t),nParticle);
x2_fw = zeros(length(t),nParticle);
v2_fw = zeros(length(t),nParticle);
x3_fw = zeros(length(t),nParticle);
v3_fw = zeros(length(t),nParticle);

% x10 = 0.5;
% v10 = 0;
% x20 = 0.5;
% v20 = 0;
% x30 = 0.5;
% v30 = 0;
x1_fw(1,:) = x1_init;
v1_fw(1,:) = v1_init;
x2_fw(1,:) = x2_init;
v2_fw(1,:) = v2_init;
x3_fw(1,:) = x3_init;
v3_fw(1,:) = v3_init;


%q = 1;
%m = 1;

f1 = [x1_fw(1,:) v1_fw(1,:)];
f2 = [x2_fw(1,:) v2_fw(1,:)];
f3 = [x3_fw(1,:) v3_fw(1,:)];

D_T = [t(1) t(end)];

% Here we calculate the forward system in the old way with ode45, to
% compare it to the VV result.
[t_fw,x_v_fwold] = forwardSystem_3D(interp_E_x, interp_E_y, interp_E_z,...
        x_v_init, D_T, nParticle, q, m);

% This is the implementation of the VV method.
for i = 2:length(t)
    
    a1_1 = (q/m)*interp_E_x(x1_fw(i-1,:), x2_fw(i-1,:), x2_fw(i-1,:));
    a2_1 = (q/m)*interp_E_y(x1_fw(i-1,:), x2_fw(i-1,:), x2_fw(i-1,:));   
    a3_1 = (q/m)*interp_E_z(x1_fw(i-1,:), x2_fw(i-1,:), x2_fw(i-1,:));
    
    
    x1_fw(i,:) = x1_fw(i-1,:) + v1_fw(i-1,:)*delta_t + 0.5*a1_1*delta_t.^2;
    x2_fw(i,:) = x2_fw(i-1,:) + v2_fw(i-1,:)*delta_t + 0.5*a2_1*delta_t.^2;
    x3_fw(i,:) = x3_fw(i-1,:) + v3_fw(i-1,:)*delta_t + 0.5*a3_1*delta_t.^2;

    
    
    a1_2 =  (q/m)*interp_E_x(x1_fw(i,:), x2_fw(i,:), x3_fw(i,:));
    a2_2 =  (q/m)*interp_E_y(x1_fw(i,:), x2_fw(i,:), x3_fw(i,:));
    a3_2 =  (q/m)*interp_E_z(x1_fw(i,:), x2_fw(i,:), x3_fw(i,:));

    
    v1_fw(i,:) = v1_fw(i-1,:) + 0.5*(a1_1 + a1_2)*delta_t;
    v2_fw(i,:) = v2_fw(i-1,:) + 0.5*(a2_1 + a2_2)*delta_t;
    v3_fw(i,:) = v3_fw(i-1,:) + 0.5*(a3_1 + a3_2)*delta_t;

    
    f1 = [f1 0.5*a1_1*delta_t.^2 0.5*(a1_1 + a1_2)*delta_t];
    f2 = [f2 0.5*a2_1*delta_t.^2 0.5*(a2_1 + a2_2)*delta_t];
    f3 = [f3 0.5*a3_1*delta_t.^2 0.5*(a3_1 + a3_2)*delta_t];
    
    
    
    
end
% 
% figure()
% hold on
% plot(t, x1_fw)
% plot(t, v1_fw)
% figure()
% hold on
% plot(t, x2_fw)
% plot(t, v2_fw)
% figure()
% hold on
% plot(t, x3_fw)
% plot(t, v3_fw)


% Now as Matrix
%t = 1:1:3;
%delta_t = 5;

m = length(t)*2*nParticle;
m = m*3;
n = m;

i1_1 = 1:1:(length(t)*2*nParticle);
j1_1 = 1:1:(length(t)*2*nParticle);
v1_1 = ones(1,length(t)*2*nParticle);

i1_2 = (2*nParticle+1):1:((length(t)*2)*nParticle);
j1_2 = 1:1:((length(t)*2-2)*nParticle);
v1_2 = -1*ones(1,(length(t)*2-2)*nParticle);

i1_3 = [];
j1_3 = [];
for ii = 1:1:(length(t)-1)
    
i1_3_a = ((2*nParticle+1)+2*nParticle*(ii-1)):1:((3*nParticle)+2*nParticle*(ii-1));
j1_3_a = (((2*(nParticle)-1)+2*nParticle*(ii-1))-(nParticle-2)):1:((3*nParticle-2)+2*nParticle*(ii-1)-(nParticle-2));

i1_3 = [i1_3 i1_3_a];
j1_3 = [j1_3 j1_3_a];
end
v1_3 = -ones(1,(length(t)-1)*nParticle)*delta_t;

i2_1 = i1_1 + length(t)*2*nParticle;
j2_1 = j1_1 + length(t)*2*nParticle;
v2_1 = v1_1;

i2_2 = i1_2 + length(t)*2*nParticle;
j2_2 = j1_2 + length(t)*2*nParticle;
v2_2 = v1_2;

i2_3 = i1_3 + length(t)*2*nParticle;
j2_3 = j1_3 + length(t)*2*nParticle;
v2_3 = v1_3;

i3_1 = i2_1 + length(t)*2*nParticle;
j3_1 = j2_1 + length(t)*2*nParticle;
v3_1 = v1_1;

i3_2 = i2_2 + length(t)*2*nParticle;
j3_2 = j2_2 + length(t)*2*nParticle;
v3_2 = v1_2;

i3_3 = i2_3 + length(t)*2*nParticle;
j3_3 = j2_3 + length(t)*2*nParticle;
v3_3 = v2_3;

% i1 = [i1_1 i1_2 i1_3];
% j1 = [j1_1 j1_2 j1_3];
% v1 = [v1_1 v1_2 v1_3];
% S1 = sparse(i1,j1,v1,m/3,n/3);
% 
% 
% i2 = [i2_1 i2_2 i2_3];
% j2 = [j2_1 j2_2 j2_3];
% v2 = [v2_1 v2_2 v2_3];
% S2 = sparse(i2,j2,v2,m/3*2,n/3*2);
% 
% 
% i3 = [i3_1 i3_2 i3_3];
% j3 = [j3_1 j3_2 j3_3];
% v3 = [v3_1 v3_2 v3_3];
% S3 = sparse(i3,j3,v3,m,n);

i = [i1_1 i1_2 i1_3 i2_1 i2_2 i2_3 i3_1 i3_2 i3_3];
j = [j1_1 j1_2 j1_3 j2_1 j2_2 j2_3 j3_1 j3_2 j3_3];
v = [v1_1 v1_2 v1_3 v2_1 v2_2 v2_3 v3_1 v3_2 v3_3];


f = [f1 f2 f3];
S = sparse(i,j,v,m,n);

x_vec = S\f';

x_vec_split = reshape(x_vec, [2*length(t)*nParticle, 3]);

x_vec1 = x_vec_split(:,1);
x_vec2 = x_vec_split(:,2);
x_vec3 = x_vec_split(:,3);

x_fw = zeros(length(t),nParticle*2*3);
xx1 = reshape(x_vec1, [nParticle, 2*length(t)]);
xx2 = reshape(x_vec2, [nParticle, 2*length(t)]);
xx3 = reshape(x_vec3, [nParticle, 2*length(t)]);

for ii = 1:1:length(t)
    
    x_fw(ii,id_x1) = xx1(:,ii*2-1);
    x_fw(ii,id_v1) = xx1(:,ii*2);
    
    x_fw(ii,id_x2) = xx2(:,ii*2-1);
    x_fw(ii,id_v2) = xx2(:,ii*2);
    
    x_fw(ii,id_x3) = xx3(:,ii*2-1);
    x_fw(ii,id_v3) = xx3(:,ii*2);
    
end
% 
% x1_2 = x_vec(1:2:(2*length(t)-1));
% v1_2 = x_vec(2:2:(2*length(t)));
% x2_2 = x_vec( (2*length(t) + 1):2:(4*length(t)-1));
% v2_2 = x_vec( (2*length(t) + 2):2:(4*length(t)));
% x3_2 = x_vec((4*length(t) + 1):2:(6*length(t)-1));
% v3_2 = x_vec((4*length(t) +2):2:(6*length(t)));
 
% figure()
% hold on 
% plot(t, x1_2)
% plot(t, v1_2)
% figure()
% hold on 
% plot(t, x2_2)
% plot(t, v2_2)
% figure()
% hold on 
% plot(t, x3_2)
% plot(t, v3_2)

end


