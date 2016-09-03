function [result, val_vec] = multipleRestriction3D(i_x , i_y, i_z, w000, w100, w010, w110, w001, w101, w011, w111, val, Nx, Ny, Nz, nParticle)
% u = 1;
result = zeros(Nx, Ny, Nz);
%clear val_vec w_vec num_vec ii_vec
assert(all(all(w000 >= -1e-6)))
assert(all(all(w100 >= -1e-6)))
assert(all(all(w010 >= -1e-6)))
assert(all(all(w110 >= -1e-6)))
assert(all(all(w001 >= -1e-6)))
assert(all(all(w101 >= -1e-6)))
assert(all(all(w010 >= -1e-6)))
assert(all(all(w111 >= -1e-6)))
% 
% if nParticle == 1
%     
%     w000 = w000';
%     w001 = w001';
%     w010 = w010';
%     w011 = w011';
%     w100 = w100';
%     w101 = w101';
%     w110 = w110';
%     w111 = w111';
% end
%assert(all(w000+w100+w010+w110+w001+w101+w010+w111))
for jj = 1:nParticle

for ii = 1:size(i_x,1)
    
    val_ii = val(ii,jj);
    %{
%     % Goal: 7,7,3
%     if( i_x(ii) == 7)
%         if (i_y(ii) == 7)
%             if(i_z(ii) == 3)
%                 disp('in if 7,7,3')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w000(ii);
%                 num_vec(u) = 773;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
%     
%     if( i_x(ii) == 6)
%         if (i_y(ii) == 7)
%             if(i_z(ii) == 3)
%                 disp('in if 6,7,3')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w100(ii);
%                 num_vec(u) = 673;
%                 u = u+1;
%                 ii_vec(u) = ii;                
%             end
%         end
%     end
%     if( i_x(ii) == 7)
%         if (i_y(ii) == 6)
%             if(i_z(ii) == 3)
%                 disp('in if 7,6,3')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w010(ii);
%                 num_vec(u) = 763;
%                 u = u+1;
%                 ii_vec(u) = ii;                 
%             end
%         end
%     end
%     if( i_x(ii) == 7)
%         if (i_y(ii) == 7)
%             if(i_z(ii) == 2)
%                 disp('in if 7,7,2')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w001(ii);
%                 num_vec(u) = 772;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
% % Goal: 7,7,3
%     
%     if( i_x(ii) == 7)
%         if (i_y(ii) == 6)
%             if(i_z(ii) == 2)
%                 disp('in if 7,6,2')                
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w011(ii);
%                 num_vec(u) = 762;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
% 
%     if( i_x(ii) == 6)
%         if (i_y(ii) == 7)
%             if(i_z(ii) == 2)
%                 disp('in if 6,7,2')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w101(ii);
%                 num_vec(u) = 672;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
%     
%     if( i_x(ii) == 6)
%         if (i_y(ii) == 6)
%             if(i_z(ii) == 3)
%                 disp('in if 6,6,3')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w110(ii);
%                 num_vec(u) = 663;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
%     
%     if( i_x(ii) == 6)
%         if (i_y(ii) == 6)
%             if(i_z(ii) == 2)
%                 disp('in if 6,6,2')
%                 val_vec(u) = val_ii;
%                 w_vec(u) = w111(ii);
%                 num_vec(u) = 662;
%                 u = u+1;
%                 ii_vec(u) = ii;
%             end
%         end
%     end
%}    
    
        result( i_x(ii,jj)      , i_y(ii,jj)       , i_z(ii,jj))        = result( i_x(ii,jj)       , i_y(ii,jj)       , i_z(ii,jj)) + w000(jj,ii)*val_ii;    
    if(i_x(ii,jj) < Nx)
        result( i_x(ii,jj) + 1  , i_y(ii,jj)       , i_z(ii,jj))        = result( i_x(ii,jj) + 1   , i_y(ii,jj)       , i_z(ii,jj)) + w100(jj,ii)*val_ii;    end    
    if(i_y(ii,jj) < Ny) 
        result( i_x(ii,jj)      , i_y(ii,jj) + 1   , i_z(ii,jj))        = result( i_x(ii,jj)       , i_y(ii,jj) + 1   , i_z(ii,jj)) + w010(jj,ii)*val_ii;    end
    if((i_x(ii,jj) < Nx) && (i_y(ii,jj) < Ny))
        result( i_x(ii,jj) + 1  , i_y(ii,jj) + 1   , i_z(ii,jj))        = result( i_x(ii,jj) + 1   , i_y(ii,jj) + 1   , i_z(ii,jj)) + w110(jj,ii)*val_ii;    end
    
    if(i_z(ii,jj) < Nz)
        result( i_x(ii,jj)      , i_y(ii,jj)       , i_z(ii,jj) + 1)    = result( i_x(ii,jj)      , i_y(ii,jj)       , i_z(ii,jj) + 1) + w001(jj,ii)*val_ii; end
    if((i_x(ii,jj) < Nx) && (i_z(ii,jj) < Nz))
        result( i_x(ii,jj) + 1  , i_y(ii,jj)       , i_z(ii,jj) + 1)    = result( i_x(ii,jj) + 1  , i_y(ii,jj)       , i_z(ii,jj) + 1) + w101(jj,ii)*val_ii; end
    if((i_y(ii,jj) < Ny) && (i_z(ii,jj) < Nz))
        result( i_x(ii,jj)      , i_y(ii,jj) + 1   , i_z(ii,jj) + 1)    = result( i_x(ii,jj)      , i_y(ii,jj) + 1   , i_z(ii,jj) + 1) + w011(jj,ii)*val_ii; end
    if((i_x(ii,jj) < Nx) && (i_y(ii,jj) < Ny) && (i_z(ii,jj) < Nz))
        result( i_x(ii,jj) + 1  , i_y(ii,jj) + 1   , i_z(ii,jj) + 1)    = result( i_x(ii,jj) + 1  , i_y(ii,jj) + 1   , i_z(ii,jj) + 1) + w111(jj,ii)*val_ii; end

    
end

end
end