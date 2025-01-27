%%
%copied from 
% https://github.com/nsteinme/nickBox/blob/master/matlab/affineMat.m

classdef affineMat
    % e.g. from https://uk.mathworks.com/discovery/affine-transformation.html
    % or https://en.wikipedia.org/wiki/Transformation_matrix
    % (note these two have different convention on rotations)
    properties
        
    end
    
    methods (Static)
        function T = trans(tx, ty)
            
            T = [1 0 tx; ...
                0 1 ty; ...
                0 0 1];
        end
        
        
        function T = rot(theta)
            
            T = [cos(theta) -sin(theta) 0;...
                sin(theta) cos(theta) 0;...
                0 0 1];
        end
        
        function T = scale(sx, sy)
            
            T = [sx 0 0; ...
                0 sy 0;...
                0 0 1];
        end
        
        function T = shear(shx, shy)
            
            T = [1 shy 0; ...
                shx 1 0;...
                0 0 1];
        end
    end
    
end