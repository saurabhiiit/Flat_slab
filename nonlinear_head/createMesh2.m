function [ nodal_connect, nodal_coordinate, faces, slab_mesh_meta_data , col_mesh_meta_data, col_bot_node, slab_side_nodes,slab_top_nodes  ] = createMesh2(slab_dimension,col_dimension,slab_divisions,col_divisions,drop_divisions,col_top_postion,drop_dimension,drop_top_postion )
%CREATEMESH Summary of this function goes here
%   Detailed explanation goes here


% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
slab_mesh_size = 50 * floor((slab_dimension./slab_divisions)/50); 
% disp(slab_mesh_size)

% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
col_mesh_size = 50 * floor((col_dimension./col_divisions)/50); 

% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
drop_mesh_size = 100 * floor((drop_dimension./drop_divisions)/100); 


%% Columns
col_x_coord1 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(1,1) - (col_dimension(1)/2)));
% disp(col_x_coord1)
col_y_coord1 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(1,2) - (col_dimension(2)/2));
col_z_coord1 = fliplr(helper(col_mesh_size(3), col_dimension(3), 0));
col_nodal_coordinate1 = combvec(col_x_coord1, col_y_coord1, col_z_coord1).';
% disp(col_nodal_coordinate1)
% disp(fliplr(col_z_coord1));
col_x_coord2 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(2,1) - (col_dimension(1)/2)));
col_y_coord2 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(2,2) - (col_dimension(2)/2));
col_nodal_coordinate2 = combvec(col_x_coord2, col_y_coord2, col_z_coord1).';

col_x_coord3 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(3,1) - (col_dimension(1)/2)));
col_y_coord3 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(3,2) - (col_dimension(2)/2));
col_nodal_coordinate3 = combvec(col_x_coord3, col_y_coord3, col_z_coord1).';

col_x_coord4 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(4,1) - (col_dimension(1)/2)));
col_y_coord4 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(4,2) - (col_dimension(2)/2));
col_nodal_coordinate4 = combvec(col_x_coord4, col_y_coord4, col_z_coord1).';

% disp(combvec(col_x_coord(:,1), col_y_coord(:,1), col_z_coord_reverse(:,1)).')
% disp(col_nodal_coordinate1)
%% Drop Panel

% drop_x_coord1 = helper(drop_mesh_size(1), drop_dimension(1),(drop_top_postion(1,1) - (drop_dimension(1)/2)));
drop_x_coord1  = [(drop_top_postion(1,1) - (drop_dimension(1)/2)) (drop_top_postion(1,1) - (drop_dimension(1)/2))+drop_dimension(1)];
% disp(drop_x_coord1)
drop_x_coord1 = helper2(col_x_coord1, drop_x_coord1);
 % disp(drop_x_coord1)
% drop_y_coord1 = helper(drop_mesh_size(2), drop_dimension(2),drop_top_postion(1,2) - (drop_dimension(2)/2));
drop_y_coord1  = [(drop_top_postion(1,2) - (drop_dimension(2)/2)) (drop_top_postion(1,2) - (drop_dimension(2)/2))+drop_dimension(2)];
drop_y_coord1 = helper2(col_y_coord1, drop_y_coord1);

drop_z_coord1 = fliplr(helper(drop_mesh_size(3), drop_dimension(3), col_top_postion(1,3)));
drop_nodal_coordinate1 = combvec(drop_x_coord1, drop_y_coord1, drop_z_coord1).';


% drop_x_coord2 = helper(drop_mesh_size(1), drop_dimension(1),(drop_top_postion(2,1) - (drop_dimension(1)/2)));
% disp(drop_x_coord1)
drop_x_coord2  = [(drop_top_postion(2,1) - (drop_dimension(1)/2)) (drop_top_postion(2,1) - (drop_dimension(1)/2))+drop_dimension(1)];

drop_x_coord2 = helper2(col_x_coord2, drop_x_coord2);
% drop_y_coord2 = helper(drop_mesh_size(2), drop_dimension(2),drop_top_postion(2,2) - (drop_dimension(2)/2));
drop_y_coord2  = [(drop_top_postion(2,2) - (drop_dimension(2)/2)) (drop_top_postion(2,2) - (drop_dimension(2)/2))+drop_dimension(2)];

drop_y_coord2 = helper2(col_y_coord2, drop_y_coord2);
drop_nodal_coordinate2 = combvec(drop_x_coord2, drop_y_coord2, drop_z_coord1).';

% drop_x_coord3 = helper(drop_mesh_size(1), drop_dimension(1),(drop_top_postion(3,1) - (drop_dimension(1)/2)));
drop_x_coord3  = [(drop_top_postion(3,1) - (drop_dimension(1)/2)) (drop_top_postion(3,1) - (drop_dimension(1)/2))+drop_dimension(1)];
drop_x_coord3 = helper2(col_x_coord3, drop_x_coord3);
% drop_y_coord3 = helper(drop_mesh_size(2), drop_dimension(2),drop_top_postion(3,2) - (drop_dimension(2)/2));
drop_y_coord3  = [(drop_top_postion(3,2) - (drop_dimension(2)/2)) (drop_top_postion(3,2) - (drop_dimension(2)/2))+drop_dimension(2)];

drop_y_coord3 = helper2(col_y_coord3, drop_y_coord3);
drop_nodal_coordinate3 = combvec(drop_x_coord3, drop_y_coord3, drop_z_coord1).';

% drop_x_coord4 = helper(drop_mesh_size(1), drop_dimension(1),(drop_top_postion(4,1) - (drop_dimension(1)/2)));
drop_x_coord4  = [(drop_top_postion(4,1) - (drop_dimension(1)/2)) (drop_top_postion(4,1) - (drop_dimension(1)/2))+drop_dimension(1)];
drop_x_coord4 = helper2(col_x_coord4, drop_x_coord4);
% drop_y_coord4 = helper(drop_mesh_size(2), drop_dimension(2),drop_top_postion(4,2) - (drop_dimension(2)/2));
drop_y_coord4  = [(drop_top_postion(4,2) - (drop_dimension(2)/2)) (drop_top_postion(4,2) - (drop_dimension(2)/2))+drop_dimension(2)];
drop_y_coord4 = helper2(col_y_coord4, drop_y_coord4);
drop_nodal_coordinate4 = combvec(drop_x_coord4, drop_y_coord4, drop_z_coord1).';

% disp(drop_nodal_coordinate1)

%% For slab

% % For slab
slab_x_coord = helper(slab_mesh_size(1), slab_dimension(1),0);
slab_x_coord = helper2(drop_x_coord1, slab_x_coord);
slab_x_coord = helper2(drop_x_coord2, slab_x_coord);

slab_y_coord = helper(slab_mesh_size(2), slab_dimension(2),0);
slab_y_coord = helper2(drop_y_coord1, slab_y_coord);
slab_y_coord = helper2(drop_y_coord3, slab_y_coord);

slab_z_coord = helper(slab_mesh_size(3), slab_dimension(3), drop_top_postion(1,3));
slab_nodal_coordinate = combvec(slab_x_coord, slab_y_coord, slab_z_coord).';
% [slab_x_coord , coordx] = slab_helper_x(slab_mesh_size(1), slab_dimension(1), 0, col_mesh_size , col_dimension,col_top_postion);
% [slab_y_coord , coordy] = slab_helper_y(slab_mesh_size(2), slab_dimension(2), 0, col_mesh_size , col_dimension,col_top_postion);
% slab_z_coord = helper(slab_mesh_size(3), slab_dimension(3), col_top_postion(1,3));	
% slab_nodal_coordinate = combvec(slab_x_coord, slab_y_coord, slab_z_coord).';

%Calculating the attached coordinate of col with slab
% disp(length(coordx(4)));

drop_top_coordinate = combvec([drop_x_coord1 drop_x_coord2], [drop_y_coord1 drop_y_coord3] ,drop_dimension(3)+col_dimension(3)).';
% drop_bot_coordinate = combvec([drop_x_coord1 drop_x_coord2], [drop_y_coord1 drop_y_coord3] ,col_dimension(3)).';
drop_bot_coordinate = [(combvec(drop_x_coord1, drop_y_coord1 ,col_dimension(3)).'); (combvec(drop_x_coord2,drop_y_coord2 ,col_dimension(3)).') ; (combvec(drop_x_coord3, drop_y_coord3 ,col_dimension(3)).');(combvec(drop_x_coord4, drop_y_coord4 ,col_dimension(3)).') ];
% disp(drop_bot_coordinate)

col_top_coordinate = combvec([col_x_coord1 col_x_coord2], [col_y_coord1 col_y_coord3] ,col_dimension(3)).';
 % disp(col_top_coordinate)

% disp(col_nodal_coordinate(:,3,1))


% For creating the nodal connectivity matrix for the slab

% Calculating the number of element in the slab
slab_no_elements = (length(slab_x_coord) - 1) * (length(slab_y_coord) - 1) * (length(slab_z_coord) - 1);

slab_mesh_meta_data = [length(slab_x_coord) - 1, length(slab_y_coord) - 1, length(slab_z_coord) - 1];
% disp(slab_no_elements)
% disp(slab_mesh_meta_data)
% slab_no_elements = 5*4*3;

% slab_mesh_meta_data = [5, 4, 3];
slab_side_nodes =[];
slab_top_nodes = [];
for ii = 1:slab_no_elements
    element_no = ii;  
    layer = floor((element_no-1)/(slab_mesh_meta_data(1)*slab_mesh_meta_data(2)));
    % disp(layer)
    temp = mod(element_no, slab_mesh_meta_data(1)*slab_mesh_meta_data(2));
      % disp(temp);
    if temp
        element_no = temp;
        % disp(temp)
    else
        element_no = slab_mesh_meta_data(1)*slab_mesh_meta_data(2);
   		% disp(element_no)
    end
    row = floor((element_no-1)/slab_mesh_meta_data(1));
     % disp(row)
    temp2 = mod(element_no, slab_mesh_meta_data(1));
    % disp(temp2)
    if temp2
        index = temp2;
    else
        index = slab_mesh_meta_data(1);
    end
    % disp(index)
    % disp(index)
    %% Mapping of local nodes to global nodes
    % First node of element. It is mapped to corresponding x node in global.
    node_one = layer*(slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1) + row*(slab_mesh_meta_data(1)+1) + index;
     % disp(node_one);
    node_two = node_one + (slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1);
   %    if layer==0
   %     disp(node_one)
   %     disp(node_two)
   % end
    mapping = [node_one, node_one + 1, node_one + slab_mesh_meta_data(1) + 2,  node_one + slab_mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + slab_mesh_meta_data(1) + 2, node_two + slab_mesh_meta_data(1) + 1];
    % disp(mapping)
    face = [mapping(1) mapping(2) mapping(3) mapping(4);
            mapping(5) mapping(6) mapping(7) mapping(8);
            mapping(1) mapping(2) mapping(6) mapping(5);
            mapping(4) mapping(3) mapping(7) mapping(8);
            mapping(2) mapping(3) mapping(7) mapping(6);
            mapping(1) mapping(4) mapping(8) mapping(5)];
    slab_nodal_connect(ii, :) = mapping;
    slab_faces(6*(ii-1)+1:6*ii, :) = face;
    %for calculating nodes of the side face of the slab
    if index==1
        slab_side_nodes=[slab_side_nodes ii];
    end
    if layer == slab_mesh_meta_data(3)-1
         slab_top_nodes =[slab_top_nodes ii];
    end

end
% slab_side_nodes = unique(slab_side_nodes);

% disp(slab_nodal_connect)
% disp((slab_nodal_connect))
 % disp(slab_faces)


% For creating nodal connectivity matrix adding the column nodal connect and coordinates

% Number of element in the column

drop_top_index = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, drop_top_coordinate);
% disp(drop_top_index)


% %Caluclate the index of the column and slab common connection node
% for ij = 1 : length(slab_x_coord)*length(slab_y_coord)
% 	for ik =1 : length(col_top_coordinate)
% 		if isequal(col_top_coordinate(ik,:),slab_nodal_coordinate(ij,:)) == 1
% 			  % disp(ij);
% 			 % disp(col_top_coordinate(ik,:))
% 			col_top_index(ik) = ij;
% 		end
% 	end
% end

drop_no_elements =  (length(drop_x_coord1) - 1) * (length(drop_y_coord1) - 1) * (length(drop_z_coord1) - 1);
drop_mesh_meta_data = [length(drop_x_coord1) - 1, length(drop_y_coord1) - 1, length(drop_z_coord1) - 1];

% disp(drop_no_elements)
% disp(drop_mesh_meta_data)

[drop_nodal_connect, drop_faces, drop_bot_node] = helper4(slab_nodal_connect, drop_no_elements, drop_mesh_meta_data , drop_top_index);
 % disp(drop_bot_node);
 % disp(length(drop_nodal_connect(:,:,3)));
for ij = 1 : length(drop_bot_coordinate)
    for ik =1 : length(col_top_coordinate)
        if isequal(col_top_coordinate(ik,:),drop_bot_coordinate(ij,:)) == 1
                % disp(ij);
             % disp(col_top_coordinate(ik,:))
            % disp(drop_bot_node(ij));
            col_top_index(ik) = drop_bot_node(ij);
             break;
        end
    end
end
 % disp(drop_bot_coordinate)
  % disp(col_top_index)
 % disp(drop_bot_node)
col_no_elements =  (length(col_x_coord1) - 1) * (length(col_y_coord1) - 1) * (length(col_z_coord1) - 1);
col_mesh_meta_data = [length(col_x_coord1) - 1, length(col_y_coord1) - 1, length(col_z_coord1) - 1];
[col_nodal_connect, col_faces, col_bot_node] = helper4(drop_nodal_connect(:,:,4), col_no_elements, col_mesh_meta_data, col_top_index);


initial_cord_col = (col_mesh_meta_data(1) + 1)*(col_mesh_meta_data(2) + 1) + 1; % This is done so that coordinate don't get repeated from slab to column
initial_cord_drop = (drop_mesh_meta_data(1) + 1)*(drop_mesh_meta_data(2) + 1) + 1;
% disp(initial_cord_col)
% Combining the nodal coordinate to the one
% nodal_coordinate = [slab_nodal_coordinate ;col_nodal_coordinate1(initial_cord_col:end,:); col_nodal_coordinate2(initial_cord_col:end,:); col_nodal_coordinate3(initial_cord_col:end,:); col_nodal_coordinate4(initial_cord_col:end,:) ];
 nodal_coordinate = [slab_nodal_coordinate ;drop_nodal_coordinate1(initial_cord_drop:end,:); drop_nodal_coordinate2(initial_cord_drop:end,:); drop_nodal_coordinate3(initial_cord_drop:end,:); drop_nodal_coordinate4(initial_cord_drop:end,:); col_nodal_coordinate1(initial_cord_col:end,:); col_nodal_coordinate2(initial_cord_col:end,:); col_nodal_coordinate3(initial_cord_col:end,:); col_nodal_coordinate4(initial_cord_col:end,:)];

    % if layer == col_mesh_meta_data(3)-1
    %              col_bot_node =[col_bot_node node_two, node_two + 1,  node_two + col_mesh_meta_data(1) + 2, node_two + col_mesh_meta_data(1) + 1];
    %         end
    %         col_bot_node = unique(col_bot_node);
faces = [slab_faces ;drop_faces(:,:,1); drop_faces(:,:,2) ;drop_faces(:,:,3); drop_faces(:,:,4); col_faces(:,:,1); col_faces(:,:,2) ;col_faces(:,:,3); col_faces(:,:,4)];
nodal_connect = [slab_nodal_connect ;drop_nodal_connect(:,:,1); drop_nodal_connect(:,:,2); drop_nodal_connect(:,:,3); drop_nodal_connect(:,:,4);col_nodal_connect(:,:,1); col_nodal_connect(:,:,2); col_nodal_connect(:,:,3); col_nodal_connect(:,:,4)];

% faces = [slab_faces ;col_faces(:,:,1); col_faces(:,:,2) ;col_faces(:,:,3); col_faces(:,:,4)];
% nodal_connect = [slab_nodal_connect ;col_nodal_connect(:,:,1); col_nodal_connect(:,:,2); col_nodal_connect(:,:,3); col_nodal_connect(:,:,4)];
 % patch('Faces',slab_faces,'Vertices',slab_nodal_coordinate,'FaceColor','blue')
 % patch('Faces',col_faces(:,:,1),'Vertices',col_nodal_coordinate1,'FaceColor','blue')
 % patch('Faces',col_faces(:,:,2),'Vertices',col_nodal_coordinate2,'FaceColor','blue')
 % patch('Faces',col_faces(:,:,3),'Vertices',col_nodal_coordinate3,'FaceColor','blue')
 % patch('Faces',col_faces(:,:,4),'Vertices',col_nodal_coordinate4,'FaceColor','blue')


% slab_side_nodes
end
function [coordinates] = helper(mesh_size, dimension,coordinates)
    flag = 0;
    while(flag == 0)
        coord_val = min(coordinates(end) + mesh_size , coordinates(1) + dimension);
        coordinates(end + 1) = coord_val;
        flag = coordinates(end) >= (coordinates(1) + dimension);
    end
end

function [coordinates] = helper2(col_x_coord, drop_x_coord)
    coordinates = [];
    for i=1:length(drop_x_coord)
        if (drop_x_coord(i) <  col_x_coord(1))|| (drop_x_coord(i) >  col_x_coord(end))
            coordinates = [coordinates col_x_coord drop_x_coord(i)];
        else
            coordinates = [coordinates col_x_coord];
        end
        coordinates = unique(coordinates);
    end  
end
function [top_index] = helper3(x_coord, y_coord, nodal_coordinate, top_coordinate)
    %Caluclate the index of the column and slab common connection node
    % disp('a');
    for ij = 1 : length(x_coord)*length(y_coord)
        for ik =1 : length(top_coordinate)
            % disp(top_coordinate(ik,:));
            if isequal(top_coordinate(ik,:),nodal_coordinate(ij,:)) == 1
                % disp(ij);
             % disp(col_top_coordinate(ik,:))
                
                top_index(ik) = ij;
                
            end
        end
    end
end

function [nodal_connect, faces, bot_node] = helper4(nodal_connectp, no_elements, mesh_meta_data, top_index)
    bot_node =[];
for noc = 1:4
    % noc =2;
for ii = 1:no_elements
    element_no = ii;  
    layer = floor((element_no-1)/(mesh_meta_data(1)*mesh_meta_data(2)));
    % disp(layer)
    temp = mod(element_no, mesh_meta_data(1)*mesh_meta_data(2));
      % disp(temp);
    if temp
        element_no = temp;
        % disp(temp)
    else
        element_no = mesh_meta_data(1)*mesh_meta_data(2);
        % disp(element_no)
    end
    row = floor((element_no-1)/mesh_meta_data(1));
     % disp(row)
    temp2 = mod(element_no, mesh_meta_data(1));
    % disp(temp2)
    if temp2
        index = temp2;
    else
        index = mesh_meta_data(1);
    end
    % disp(layer);    

    if layer == 0
        if noc == 1 
        node_one = top_index(2*row*(mesh_meta_data(1)+1) + index);
        node_two = max(nodal_connectp(:)) + row*(mesh_meta_data(1)+1) + index;
         %  disp(node_one);
         % disp(node_two);
        mapping = [node_one , node_one + 1 , top_index(2*(row+1)*(mesh_meta_data(1)+1) + index +  1) , top_index(2*(row+1)*(mesh_meta_data(1)+1) + index ) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
        col1max = max(mapping(:));
        if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
        % disp(mapping);
        % face = [mapping(1) mapping(2) mapping(3) mapping(4);
     %        mapping(5) mapping(6) mapping(7) mapping(8);
     %        mapping(1) mapping(2) mapping(6) mapping(5);
     %        mapping(4) mapping(3) mapping(7) mapping(8);
     %        mapping(2) mapping(3) mapping(7) mapping(6);
     %        mapping(1) mapping(4) mapping(8) mapping(5)];
        % nodal_connect(ii, :) = mapping;
        % faces(6*(ii-1)+1:6*ii, :) = face;

        elseif noc == 2
            node_one = top_index(2*(row)*(mesh_meta_data(1)+1) + index + (mesh_meta_data(1)+1));
            node_two = col1max + row*(mesh_meta_data(1)+1) + index;
         %  disp(node_one);
         % disp(node_two);
            mapping = [node_one , node_one + 1 , top_index(2*(row+1)*(mesh_meta_data(1)+1) + index +  1 + (mesh_meta_data(1)+1)) , top_index(2*(row+1)*(mesh_meta_data(1)+1) + index + (mesh_meta_data(1)+1)) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
            col2max = max(mapping(:));
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
            % disp(mapping)
        elseif noc == 3
            node_one = top_index(2*(mesh_meta_data(2) + 1 )*(mesh_meta_data(1) + 1 ) + 2*row*(mesh_meta_data(1) + 1 ) + index);
            node_two = col2max + row*(mesh_meta_data(1)+1) + index;
            % disp(node_one)
            mapping = [node_one , node_one + 1 , top_index(2*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) +  2*(row + 1)*(mesh_meta_data(1)+1) + index + 1 ) , top_index(2*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + 2*(row + 1)*(mesh_meta_data(1)+1) + index )  , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
            col3max = max(mapping(:));
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
              % disp(mapping)
        elseif noc == 4
            node_one = top_index(2*(mesh_meta_data(2) + 1 )*(mesh_meta_data(1) + 1 ) + 2*(row)*(mesh_meta_data(1)+1) + index + (mesh_meta_data(1)+1));
            node_two = col3max + row*(mesh_meta_data(1)+1) + index;
         %  disp(node_one);
         % disp(node_two);
            mapping = [node_one , node_one + 1 , top_index(2*(mesh_meta_data(2) + 1 )*(mesh_meta_data(1) + 1 ) + 2*(row+1)*(mesh_meta_data(1)+1) + index +  1 + (mesh_meta_data(1)+1)) , top_index(2*(mesh_meta_data(2) + 1 )*(mesh_meta_data(1) + 1 ) + 2*(row+1)*(mesh_meta_data(1)+1) + index + (mesh_meta_data(1)+1)) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);

            % disp(mapping)
        end
    else
        if noc==1
            node_one = max(nodal_connectp(:)) + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
            % disp(noc);
        % disp(node_one)
        % disp(node_two)
            mapping = [node_one, node_one + 1, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            col1max = max(mapping(:));
            % disp(col1max);
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
        elseif noc == 2
            node_one = col1max + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
         % disp(node_one)
         % disp(node_two)
            mapping = [node_one, node_one + 1, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
        % disp(mapping)
            col2max = max(mapping(:));
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
        elseif noc == 3
            node_one = col2max + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
         % disp(node_one)
         % disp(node_two)
            mapping = [node_one, node_one + 1, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
        % disp(mapping)
            col3max = max(mapping(:));
           if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
        elseif noc == 4
            node_one = col3max + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
         % disp(node_one)
         % disp(node_two)
            mapping = [node_one, node_one + 1, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
        % disp(mapping)
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
        end

        % disp(mapping)
        
    end
        face = [mapping(1) mapping(2) mapping(3) mapping(4);
            mapping(5) mapping(6) mapping(7) mapping(8);
            mapping(1) mapping(2) mapping(6) mapping(5);
            mapping(4) mapping(3) mapping(7) mapping(8);
            mapping(2) mapping(3) mapping(7) mapping(6);
            mapping(1) mapping(4) mapping(8) mapping(5)];
        nodal_connect(ii, :,noc) = mapping;
        faces(6*(ii-1)+1:6*ii, :,noc) = face;
end
end
end
