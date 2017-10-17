function [ nodal_connect, nodal_coordinate, faces, slab_mesh_meta_data , col_mesh_meta_data, col_bot_node, slab_side_nodes,slab_top_nodes,punching_element  ] = createMesh(slab_dimension,col_dimension,slab_divisions,col_divisions,col_top_postion )
%CREATEMESH Summary of this function goes here
%   Detailed explanation goes here


% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
slab_mesh_size = 50 * floor((slab_dimension./slab_divisions)/50); 
% disp(slab_mesh_size)

% Gives x_dim, y_dim and z_dim i.e. the mesh size in all the directions.
col_mesh_size = 50 * floor((col_dimension./col_divisions)/50); 



%%% column mesh division

col_mesh_size = 50 * floor((col_dimension./col_divisions)/50); 

% col 1
col_x_coord1 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(1,1) - (col_dimension(1)/2)));
col_y_coord1 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(1,2) - (col_dimension(2)/2));  
col_z_coord1 = fliplr(helper(col_mesh_size(3), col_dimension(3), 0));

% col2
col_x_coord2 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(2,1) - (col_dimension(1)/2)));
col_y_coord2 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(2,2) - (col_dimension(2)/2));

col_x_coord3 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(3,1) - (col_dimension(1)/2)));
col_y_coord3 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(3,2) - (col_dimension(2)/2));

col_x_coord4 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(4,1) - (col_dimension(1)/2)));
col_y_coord4 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(4,2) - (col_dimension(2)/2));

%%% Combining the divisons in the column to column coord to form the combined matrix
col_nodal_coordinate1 = combvec(col_x_coord1, col_y_coord1, col_z_coord1).';
col_nodal_coordinate2 = combvec(col_x_coord2, col_y_coord2, col_z_coord1).';
col_nodal_coordinate3 = combvec(col_x_coord3, col_y_coord3, col_z_coord1).';
col_nodal_coordinate4 = combvec(col_x_coord4, col_y_coord4, col_z_coord1).';


%%% Combining the divisons in the slab coord to form the combined matrix

%%% Slab mesh division .
slab_mesh_size = 50 * floor((slab_dimension./slab_divisions)/50); 

slab_x_coord = helper(slab_mesh_size(1), slab_dimension(1),0);
slab_y_coord = helper(slab_mesh_size(2), slab_dimension(2),0);
% slab_z_coord = helper(slab_mesh_size(3), slab_dimension(3), col_top_postion(1,3));
slab_z_coord = [col_top_postion(1,3) slab_dimension(3)+col_top_postion(1,3)];
% slab_x_coord = [slab_x_coord 2300 ];
% slab_y_coord = [slab_y_coord 2300 ];

% slab_x_coord  = unique([slab_x_coord col_x_coord1 col_x_coord2 ]);
 slab_x_coord  = unique([slab_x_coord(slab_x_coord<col_x_coord1(1)|slab_x_coord>col_x_coord2(end)|(slab_x_coord>col_x_coord1(end)&slab_x_coord<col_x_coord2(1))) col_x_coord1 col_x_coord2 ]);

% slab_y_coord  = unique([slab_y_coord col_y_coord1 col_y_coord3 ]);
slab_y_coord  = unique([slab_y_coord(slab_y_coord<col_y_coord1(1)|slab_y_coord>col_y_coord3(end)|(slab_y_coord>col_y_coord1(end)&slab_y_coord<col_y_coord3(1))) col_y_coord1 col_y_coord3 ]);

% for ii = 1:length(slab_x_coord)
%     if slab_x_coord(ii) == 2000
%         slab_x_coord(ii) = 1250;
%     elseif slab_x_coord(ii) == 3000
%         slab_x_coord(ii) = 3750; 
%     elseif slab_x_coord(ii) == 12000
%             slab_x_coord(ii) = 11250; 
%         elseif slab_x_coord(ii) == 13000
%             slab_x_coord(ii) = 13750; 
%     end
%       if slab_y_coord(ii) == 2000
%         slab_y_coord(ii) = 1250;
%     elseif slab_y_coord(ii) == 3000
%         slab_y_coord(ii) = 3750; 
%     elseif slab_y_coord(ii) == 12000
%             slab_y_coord(ii) = 11250; 
%         elseif slab_y_coord(ii) == 13000
%             slab_y_coord(ii) = 13750; 
%       end
% end
% slab_x_coord
slab_nodal_coordinate = combvec(slab_x_coord, slab_y_coord, slab_z_coord).';
%%%

%%% Now to find out the nodal coordinate value of the slab

[slab_nodal_connect , slab_faces,  slab_mesh_meta_data, slab_side_nodes,slab_top_nodes] = helper7(slab_x_coord, slab_y_coord, slab_z_coord);

%%%

%%% Top coordinate of the column
col_top_coordinate1 = combvec(col_x_coord1 , col_y_coord1  ,col_dimension(3)).';
col_top_coordinate2 = combvec( col_x_coord2, col_y_coord2  ,col_dimension(3)).';
col_top_coordinate3 = combvec(col_x_coord3 , col_y_coord3 ,col_dimension(3)).';
col_top_coordinate4 = combvec(col_x_coord4 ,  col_y_coord4 ,col_dimension(3)).';

col_top_coordinate(:,:,1)= col_top_coordinate1;
col_top_coordinate(:,:,2)= col_top_coordinate2;
col_top_coordinate(:,:,2)= col_top_coordinate3;
col_top_coordinate(:,:,4)= col_top_coordinate4;



%%% Top index of the column
col_top_index1 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate1);
col_top_index2 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate2);
col_top_index3 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate3);
col_top_index4 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate4);

%%% Nodal connectivity of the each four column
col_bot_node=[];
[col_nodal_connect1, col_faces1, col_mesh_meta_data1, col_bot_node1, col_last_index1] = helper4( max(slab_nodal_connect(:)), col_x_coord1, col_y_coord1, col_z_coord1, col_top_index1,col_bot_node);
[col_nodal_connect2, col_faces2, col_mesh_meta_data2, col_bot_node2, col_last_index2] = helper4( col_last_index1, col_x_coord2, col_y_coord2, col_z_coord1, col_top_index2,col_bot_node);
[col_nodal_connect3, col_faces3, col_mesh_meta_data3,col_bot_node3, col_last_index3] = helper4( col_last_index2, col_x_coord3, col_y_coord3, col_z_coord1, col_top_index3,col_bot_node);
[col_nodal_connect4, col_faces4, col_mesh_meta_data4,col_bot_node4, col_last_index4] = helper4( col_last_index3, col_x_coord4, col_y_coord4, col_z_coord1, col_top_index4,col_bot_node);
col_bot_node = [col_bot_node1 col_bot_node2 col_bot_node3 col_bot_node4];
bot_node_len = [length(col_bot_node1) length(col_bot_node2) length(col_bot_node3) length(col_bot_node4)];

col_mesh_meta_data = [col_mesh_meta_data1;col_mesh_meta_data2;col_mesh_meta_data3;col_mesh_meta_data4];
%%%


%%% Combining the nodal coordinate of slab and all columns to one
nodal_coordinate = [slab_nodal_coordinate ;col_nodal_coordinate1(length(col_top_index1)+1:end,:); col_nodal_coordinate2(length(col_top_index2)+1:end,:); col_nodal_coordinate3(length(col_top_index3)+1:end,:); col_nodal_coordinate4(length(col_top_index4)+1:end,:) ];
faces = [slab_faces ; col_faces1; col_faces2 ;col_faces3; col_faces4];
nodal_connect = [slab_nodal_connect ;col_nodal_connect1; col_nodal_connect2; col_nodal_connect3; col_nodal_connect4];
%%%

%%%punching values
% [element1]= helper9(slab_nodal_connect, col_top_index1,slab_mesh_meta_data);
% [element2]= helper9(slab_nodal_connect, col_top_index2,slab_mesh_meta_data);
% [element3]= helper9(slab_nodal_connect, col_top_index3,slab_mesh_meta_data);
% [element4]= helper9(slab_nodal_connect, col_top_index4,slab_mesh_meta_data);
% punching_element = [element1;element2;element3;element4];
punching_element =1;
end
function [coordinates] = helper(mesh_size, dimension,coordinates)
    flag = 0;
    while(flag == 0)
        coord_val = min(coordinates(end) + mesh_size , coordinates(1) + dimension);
        coordinates(end + 1) = coord_val;
        flag = coordinates(end) >= (coordinates(1) + dimension);
    end
end


function [slab_nodal_connect , slab_faces, slab_mesh_meta_data,slab_side_nodes,slab_top_nodes] = helper7(slab_x_coord, slab_y_coord, slab_z_coord)
slab_no_elements = (length(slab_x_coord) - 1) * (length(slab_y_coord) - 1) * (length(slab_z_coord) - 1);

slab_mesh_meta_data = [length(slab_x_coord) - 1, length(slab_y_coord) - 1, length(slab_z_coord) - 1];
slab_side_nodes =[];
slab_top_nodes = [];
slab_nodal_connect = zeros(slab_no_elements, 8);
slab_faces = zeros(slab_no_elements*6, 4);
for ii = 1:slab_no_elements
    element_no = ii;  
    layer = floor((element_no-1)/(slab_mesh_meta_data(1)*slab_mesh_meta_data(2)));
    temp = mod(element_no, slab_mesh_meta_data(1)*slab_mesh_meta_data(2));
    if temp
        element_no = temp;
    else
        element_no = slab_mesh_meta_data(1)*slab_mesh_meta_data(2);
    end
    row = floor((element_no-1)/slab_mesh_meta_data(1));
    temp2 = mod(element_no, slab_mesh_meta_data(1));
    if temp2
        index = temp2;
    else
        index = slab_mesh_meta_data(1);
    end
    
    %%% Mapping of local nodes to global nodes
    % First node of element. It is mapped to corresponding x node in global.
    node_one = layer*(slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1) + row*(slab_mesh_meta_data(1)+1) + index;
    node_two = node_one + (slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1);
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
% slab_top_nodes = unique(slab_top_nodes);
% slab_side_nodes = unique(slab_side_nodes);
end

function [top_index] = helper3(x_coord, y_coord, slab_nodal_coordinate, col_top_coordinate)
    %Caluclate the index of the column and slab common connection node
    ik=1;
    for ij = 1 : length(x_coord)*length(y_coord)
                    if isequal(col_top_coordinate(ik,:),slab_nodal_coordinate(ij,:)) == 1
                top_index(ik) = ij;
                 if length(col_top_coordinate)>ik
                  ik = ik + 1 ;
             end
            end
    end
end

function [nodal_connect, faces,mesh_meta_data, bot_node, last_index] = helper4(starting_index,col_x_coord,col_y_coord, col_z_coord, top_index,bot_node)
no_elements =  (length(col_x_coord) - 1) * (length(col_y_coord) - 1) * (length(col_z_coord) - 1);

mesh_meta_data = [length(col_x_coord) - 1, length(col_y_coord) - 1, length(col_z_coord) - 1]; 

    % bot_node =[];
for ii = 1:no_elements
    element_no = ii;  
    layer = floor((element_no-1)/(mesh_meta_data(1)*mesh_meta_data(2)));
    temp = mod(element_no, mesh_meta_data(1)*mesh_meta_data(2));
    if temp
        element_no = temp;
    else
        element_no = mesh_meta_data(1)*mesh_meta_data(2);
    end
    row = floor((element_no-1)/mesh_meta_data(1));
    temp2 = mod(element_no, mesh_meta_data(1));
    if temp2
        index = temp2;
    else
        index = mesh_meta_data(1);
    end
    if layer == 0
        node_one = top_index(row*(mesh_meta_data(1)+1) + index);
        node_two = starting_index + row*(mesh_meta_data(1)+1) + index;
        mapping = [node_one , node_one + 1 , top_index((row+1)*(mesh_meta_data(1)+1) + index +  1) , top_index((row+1)*(mesh_meta_data(1)+1) + index ) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
        last_index = max(mapping(:));
        if layer == mesh_meta_data(3)-1
                bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
        end
            bot_node = unique(bot_node);


    else
            node_one = starting_index + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
            mapping = [node_one, node_one + 1, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            last_index = max(mapping(:));
            if layer == mesh_meta_data(3)-1
                 bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
            end
            bot_node = unique(bot_node);
    end        
        face = [mapping(1) mapping(2) mapping(3) mapping(4);
            mapping(5) mapping(6) mapping(7) mapping(8);
            mapping(1) mapping(2) mapping(6) mapping(5);
            mapping(4) mapping(3) mapping(7) mapping(8);
            mapping(2) mapping(3) mapping(7) mapping(6);
            mapping(1) mapping(4) mapping(8) mapping(5)];
        nodal_connect(ii, :) = mapping;
        faces(6*(ii-1)+1:6*ii, :) = face;
    end
end

function [element]= helper9(slab_nodal_connect, col_top_index,slab_mesh_meta_data)
for ii=1:slab_mesh_meta_data(1)*slab_mesh_meta_data(2)
    ind1 =ismember(col_top_index,slab_nodal_connect(ii,:));
    ind2 = sum(ind1(:)==1);
    if(length(col_top_index) == ind2)
    element = ii;
    break;
    end
end
end
        
        
        