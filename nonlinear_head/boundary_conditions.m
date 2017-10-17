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



 for ii =1:length(col_bot_dof)
  global_stiff(col_bot_dof(ii),:) = 0;
  global_stiff(:,col_bot_dof(ii)) = 0;
  global_stiff(col_bot_dof(ii), col_bot_dof(ii)) =1;
  load(col_bot_dof(ii)) = 0;
 end
end