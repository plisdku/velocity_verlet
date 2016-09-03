%Dual System 
function [x_d] = dualVV(nParticle, D_G, S, t)
    g = zeros(length(t)*2*nParticle*3,1);
    [id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]  = get_Index3D(nParticle);

    % reorder D_G

    %D_G_n = [D_G(1:nParticle) D_G((3*nParticle+1):(4*nParticle)) D_G((nParticle+1):(2*nParticle)) D_G((4*nParticle+1):(5*nParticle)) D_G((2*nParticle+1):(3*nParticle)) D_G((5*nParticle+1):(6*nParticle))];
    c = length(t)*2*nParticle;
    g((c-nParticle*2+1) :c) = [D_G(1:nParticle) D_G((3*nParticle+1):(4*nParticle))];
    g((2*c-nParticle*2+1) :2*c) = [D_G((nParticle+1):(2*nParticle)) D_G((4*nParticle+1):(5*nParticle))];
    g((3*c-nParticle*2+1) :3*c) = [D_G((2*nParticle+1):(3*nParticle)) D_G((5*nParticle+1):(6*nParticle))];


    %g((end-nParticle*2*3+1):end) = D_G_n;
    %S_full = full(S);
    %S_tr = sparse(S_full');

    x_v_dual = S' \ g;

    x_vec_split = reshape(x_v_dual, [2*length(t)*nParticle, 3]);

    x_vec1 = x_vec_split(:,1);
    x_vec2 = x_vec_split(:,2);
    x_vec3 = x_vec_split(:,3);

    x_d = zeros(length(t),nParticle*2*3);
    xx1 = reshape(x_vec1, [nParticle, 2*length(t)]);
    xx2 = reshape(x_vec2, [nParticle, 2*length(t)]);
    xx3 = reshape(x_vec3, [nParticle, 2*length(t)]);

    for ii = 1:1:length(t)

        x_d(ii,id_x1) = xx1(:,ii*2-1);
        x_d(ii,id_v1) = xx1(:,ii*2);

        x_d(ii,id_x2) = xx2(:,ii*2-1);
        x_d(ii,id_v2) = xx2(:,ii*2);

        x_d(ii,id_x3) = xx3(:,ii*2-1);
        x_d(ii,id_v3) = xx3(:,ii*2);
    end

end

