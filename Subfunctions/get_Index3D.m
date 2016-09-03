function [id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]= get_Index3D(nParticle)
% [id_x1, id_x2, id_x3, id_v1, id_v2, id_v3]= get_Index3D(nParticle)
id_x1 = 1:nParticle;
id_x2 = id_x1 + nParticle;
id_x3 = id_x2 + nParticle;
id_v1 = id_x3 + nParticle;
id_v2 = id_v1 + nParticle;
id_v3 = id_v2 + nParticle;

end