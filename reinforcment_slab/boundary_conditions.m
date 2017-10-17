function [global_stiff, load] =  boundary_conditions( global_stiff,  load, col_bot_dof)
%**************************************************************************
% This function applies the given boundary condition over load and 
% stiffness matrices.
%**************************************************************************
%
% Input parameters:
% condition      - Boundary condition to apply.
% global_stiff   - Original global stiffness matrix.
% mesh_meta_data - Mesh meta data consists of number division in all the
%                  directions
% load           - Load vector applied.
%
% Output:
% global_stiff   - Global Stiffness Matrix after application of boundary
%                  conditions.
% load           - Load vector after application of BC


global_stiff(:,col_bot_dof) = 0;
global_stiff(col_bot_dof,:) = 0;
global_stiff(col_bot_dof, col_bot_dof) = eye(length(col_bot_dof));

load(col_bot_dof) = 0;
 
end