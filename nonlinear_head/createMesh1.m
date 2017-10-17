function [ nodal_connect, nodal_coordinate, faces, slab_mesh_meta_data ,col_bot_node, slab_side_nodes,slab_top_nodes  ] = createMesh1(slab_dimension,col_dimension,slab_divisions,col_divisions,drop_divisions,col_top_postion,drop_dimension,drop_top_postion )

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


%% Drop Panel

new_drop_x_coord1 = [(drop_top_postion(1,1) - (drop_dimension(1)/2)) (drop_top_postion(1,1) - (drop_dimension(1)/2))+drop_dimension(1)];
new_drop_y_coord1  = [(drop_top_postion(1,2) - (drop_dimension(2)/2)) (drop_top_postion(1,2) - (drop_dimension(2)/2))+drop_dimension(2)];
new_drop_x_coord2  = [(drop_top_postion(2,1) - (drop_dimension(1)/2)) (drop_top_postion(2,1) - (drop_dimension(1)/2))+drop_dimension(1)];
new_drop_y_coord2  = [(drop_top_postion(2,2) - (drop_dimension(2)/2)) (drop_top_postion(2,2) - (drop_dimension(2)/2))+drop_dimension(2)];
new_drop_x_coord3  = [(drop_top_postion(3,1) - (drop_dimension(1)/2)) (drop_top_postion(3,1) - (drop_dimension(1)/2))+drop_dimension(1)];
new_drop_y_coord3  = [(drop_top_postion(3,2) - (drop_dimension(2)/2)) (drop_top_postion(3,2) - (drop_dimension(2)/2))+drop_dimension(2)];
new_drop_x_coord4  = [(drop_top_postion(4,1) - (drop_dimension(1)/2)) (drop_top_postion(4,1) - (drop_dimension(1)/2))+drop_dimension(1)];
new_drop_y_coord4  = [(drop_top_postion(4,2) - (drop_dimension(2)/2)) (drop_top_postion(4,2) - (drop_dimension(2)/2))+drop_dimension(2)];
new_drop_z_coord1 = [drop_top_postion(1,3) drop_top_postion(1,3)-150];
new_drop_nodal_coordinate1 = combvec(new_drop_x_coord1, new_drop_y_coord1, new_drop_z_coord1).';
new_drop_nodal_coordinate2 = combvec(new_drop_x_coord2, new_drop_y_coord2, new_drop_z_coord1).';
new_drop_nodal_coordinate3 = combvec(new_drop_x_coord3, new_drop_y_coord3, new_drop_z_coord1).';
new_drop_nodal_coordinate4 = combvec(new_drop_x_coord4, new_drop_y_coord4, new_drop_z_coord1).';


% drop_x_coord1 = helper(drop_mesh_size(1), drop_dimension(1),(drop_top_postion(1,1) - (drop_dimension(1)/2)));
drop_x_coord1  = [(drop_top_postion(1,1) - (drop_dimension(1)/2)) (drop_top_postion(1,1) - (drop_dimension(1)/2))+drop_dimension(1)];
% disp(drop_x_coord1)
drop_x_coord1 = helper2(col_x_coord1, drop_x_coord1);
 % disp(drop_x_coord1)
% drop_y_coord1 = helper(drop_mesh_size(2), drop_dimension(2),drop_top_postion(1,2) - (drop_dimension(2)/2));
drop_y_coord1  = [(drop_top_postion(1,2) - (drop_dimension(2)/2)) (drop_top_postion(1,2) - (drop_dimension(2)/2))+drop_dimension(2)];
drop_y_coord1 = helper2(col_y_coord1, drop_y_coord1);

drop_z_coord1 = [new_drop_z_coord1(end) new_drop_z_coord1(end)-50];
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


%% For slab

% % For slab
slab_x_coord = helper(slab_mesh_size(1), slab_dimension(1),0);
slab_x_coord = helper2(new_drop_x_coord1, slab_x_coord);
slab_x_coord = helper2(new_drop_x_coord2, slab_x_coord);

slab_y_coord = helper(slab_mesh_size(2), slab_dimension(2),0);
slab_y_coord = helper2(new_drop_y_coord1, slab_y_coord);
slab_y_coord = helper2(new_drop_y_coord3, slab_y_coord);

slab_z_coord = helper(slab_mesh_size(3), slab_dimension(3), drop_top_postion(1,3));
slab_nodal_coordinate = combvec(slab_x_coord, slab_y_coord, slab_z_coord).';


[slab_nodal_connect , slab_faces,  slab_mesh_meta_data, slab_side_nodes,slab_top_nodes] = helper7(slab_x_coord, slab_y_coord, slab_z_coord);


%%% Top coordinate of the new drop
new_drop_top_coordinate1 = combvec(new_drop_x_coord1 , new_drop_y_coord1  ,col_dimension(3)+drop_dimension(3)).';
new_drop_top_coordinate2 = combvec( new_drop_x_coord2, new_drop_y_coord2  ,col_dimension(3)+drop_dimension(3)).';
new_drop_top_coordinate3 = combvec(new_drop_x_coord3 , new_drop_y_coord3 ,col_dimension(3)+drop_dimension(3)).';
new_drop_top_coordinate4 = combvec(new_drop_x_coord4 ,  new_drop_y_coord4 ,col_dimension(3)+drop_dimension(3)).';

%%% Top index of the new drop
new_drop_top_index1 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, new_drop_top_coordinate1);
new_drop_top_index2 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, new_drop_top_coordinate2);
new_drop_top_index3 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, new_drop_top_coordinate3);
new_drop_top_index4 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, new_drop_top_coordinate4);

%%% Nodal connectivity of the each four column
new_drop_bot_node=[];
% max(slab_nodal_connect(:))
[new_drop_nodal_connect1, new_drop_faces1, new_drop_bot_node1, new_drop_last_index1] = helper4( max(slab_nodal_connect(:)), new_drop_x_coord1, new_drop_y_coord1, new_drop_z_coord1, new_drop_top_index1,new_drop_bot_node);
[new_drop_nodal_connect2, new_drop_faces2, new_drop_bot_node2, new_drop_last_index2] = helper4( new_drop_last_index1, new_drop_x_coord2, new_drop_y_coord2, new_drop_z_coord1, new_drop_top_index2,new_drop_bot_node);
[new_drop_nodal_connect3, new_drop_faces3, new_drop_bot_node3, new_drop_last_index3] = helper4( new_drop_last_index2, new_drop_x_coord3, new_drop_y_coord3, new_drop_z_coord1, new_drop_top_index3,new_drop_bot_node);
[new_drop_nodal_connect4, new_drop_faces4, new_drop_bot_node4, new_drop_last_index4] = helper4( new_drop_last_index3, new_drop_x_coord4, new_drop_y_coord4, new_drop_z_coord1, new_drop_top_index4,new_drop_bot_node);
new_drop_bot_node = [new_drop_bot_node1 new_drop_bot_node2 new_drop_bot_node3 new_drop_bot_node4];


drop_bot_node=[];

max1 = max(new_drop_bot_node4);
% max1 = new_drop_last_index1;
drop_top_index1 = [new_drop_bot_node1(1),max1+1,max1+2, new_drop_bot_node1(2),max1+3,max1+4,max1+5,max1+6,max1+7,max1+8,max1+9,max1+10,new_drop_bot_node1(3),max1+11,max1+12,new_drop_bot_node1(4)];
max1 = max(drop_top_index1);
 [drop_nodal_connect1, drop_faces1, drop_bot_node1, drop_last_index1] = helper5( max1, drop_x_coord1, drop_y_coord1, drop_z_coord1, drop_top_index1,drop_bot_node);
max1=drop_last_index1;
drop_top_index2 = [new_drop_bot_node2(1),max1+1,max1+2, new_drop_bot_node2(2),max1+3,max1+4,max1+5,max1+6,max1+7,max1+8,max1+9,max1+10,new_drop_bot_node2(3),max1+11,max1+12,new_drop_bot_node2(4)];
max1 = max(drop_top_index2);
[drop_nodal_connect2, drop_faces2, drop_bot_node2, drop_last_index2] = helper5( max1, drop_x_coord2, drop_y_coord2, drop_z_coord1, drop_top_index2,drop_bot_node);
max1=drop_last_index2;


drop_top_index3 = [new_drop_bot_node3(1),max1+1,max1+2, new_drop_bot_node3(2),max1+3,max1+4,max1+5,max1+6,max1+7,max1+8,max1+9,max1+10,new_drop_bot_node3(3),max1+11,max1+12,new_drop_bot_node3(4)];
max1 = max(drop_top_index3);
[drop_nodal_connect3, drop_faces3, drop_bot_node3, drop_last_index3] = helper5( max1, drop_x_coord3, drop_y_coord3, drop_z_coord1, drop_top_index3,drop_bot_node);

max1=drop_last_index3;

drop_top_index4 = [new_drop_bot_node4(1),max1+1,max1+2, new_drop_bot_node4(2),max1+3,max1+4,max1+5,max1+6,max1+7,max1+8,max1+9,max1+10,new_drop_bot_node4(3),max1+11,max1+12,new_drop_bot_node4(4)];

% drop_bot_node=[];
% [drop_nodal_connect1, drop_faces1, drop_bot_node1, drop_last_index1] = helper4( max1, drop_x_coord1, drop_y_coord1, drop_z_coord1, drop_top_index1,drop_bot_node);
% [drop_nodal_connect2, drop_faces2, drop_bot_node2, drop_last_index2] = helper4( drop_last_index1, drop_x_coord2, drop_y_coord2, drop_z_coord1, drop_top_index2,drop_bot_node);
% [drop_nodal_connect3, drop_faces3, drop_bot_node3, drop_last_index3] = helper4( drop_last_index2, drop_x_coord3, drop_y_coord3, drop_z_coord1, drop_top_index3,drop_bot_node);
max1 = max(drop_top_index4);
[drop_nodal_connect4, drop_faces4, drop_bot_node4, drop_last_index4] = helper5( max1, drop_x_coord4, drop_y_coord4, drop_z_coord1, drop_top_index4,drop_bot_node);
drop_bot_node = [drop_bot_node1 drop_bot_node2 drop_bot_node3 drop_bot_node4];


drop_bot_node_len = [length(drop_bot_node1) length(drop_bot_node2) length(drop_bot_node3) length(drop_bot_node4)];

%%%
%%% Drop bottom coordinate
drop_bot_coordinate1 = (combvec(drop_x_coord1, drop_y_coord1 ,col_dimension(3)).');
drop_bot_coordinate2 = (combvec(drop_x_coord2,drop_y_coord2 ,col_dimension(3)).') ; 
drop_bot_coordinate3 = (combvec(drop_x_coord3, drop_y_coord3 ,col_dimension(3)).');
drop_bot_coordinate4 = (combvec(drop_x_coord4, drop_y_coord4 ,col_dimension(3)).') ;


%%% Top coordinate of the column
col_top_coordinate1 = combvec(col_x_coord1 , col_y_coord1  ,col_dimension(3)).';
col_top_coordinate2 = combvec( col_x_coord2, col_y_coord2  ,col_dimension(3)).';
col_top_coordinate3 = combvec(col_x_coord3 , col_y_coord3 ,col_dimension(3)).';
col_top_coordinate4 = combvec(col_x_coord4 ,  col_y_coord4 ,col_dimension(3)).';
  
%%% Top index of the column
% col_top_index1 = (drop_last_index1-(length(drop_x_coord1)*length(drop_y_coord1)))+[helper3(drop_x_coord1, drop_y_coord1, drop_nodal_coordinate1, col_top_coordinate1)];
col_top_index1 = helper10(drop_bot_coordinate1,col_top_coordinate1,drop_bot_node1);
col_top_index2 = helper10(drop_bot_coordinate2,col_top_coordinate2,drop_bot_node2);
col_top_index3 = helper10(drop_bot_coordinate3,col_top_coordinate3,drop_bot_node3);
col_top_index4 = helper10(drop_bot_coordinate4,col_top_coordinate4,drop_bot_node4);

% col_top_index2 = (drop_last_index2-(length(drop_x_coord2)*length(drop_y_coord2)))+[helper3(drop_x_coord2, drop_y_coord2, drop_nodal_coordinate2, col_top_coordinate2)];
% col_top_index3 = (drop_last_index3-(length(drop_x_coord3)*length(drop_y_coord3)))+[helper3(drop_x_coord3, drop_y_coord3, drop_nodal_coordinate3, col_top_coordinate3)];
% col_top_index4 = (drop_last_index4-(length(drop_x_coord4)*length(drop_y_coord4)))+[helper3(drop_x_coord4, drop_y_coord4, drop_nodal_coordinate4, col_top_coordinate4)];

%%% Nodal connectivity of the each four column
col_bot_node=[];
[col_nodal_connect1, col_faces1, col_bot_node1, col_last_index1] = helper4( max(drop_nodal_connect4(:)), col_x_coord1, col_y_coord1, col_z_coord1, col_top_index1,col_bot_node);
[col_nodal_connect2, col_faces2, col_bot_node2, col_last_index2] = helper4( col_last_index1, col_x_coord2, col_y_coord2, col_z_coord1, col_top_index2,col_bot_node);
[col_nodal_connect3, col_faces3, col_bot_node3, col_last_index3] = helper4( col_last_index2, col_x_coord3, col_y_coord3, col_z_coord1, col_top_index3,col_bot_node);
[col_nodal_connect4, col_faces4, col_bot_node4, col_last_index4] = helper4( col_last_index3, col_x_coord4, col_y_coord4, col_z_coord1, col_top_index4,col_bot_node);
col_bot_node = [col_bot_node1 col_bot_node2 col_bot_node3 col_bot_node4];
% disp(col_bot_node3);
col_bot_node_len = [length(col_bot_node1) length(col_bot_node2) length(col_bot_node3) length(col_bot_node4)];


%%% Combining the nodal coordinate of slab and all drop and all columns to one

drop_nodal_index =[2 3 5 6 7 8 9 10 11 12 14 15] ;
nodal_coordinate = [slab_nodal_coordinate; 
    new_drop_nodal_coordinate1(length(new_drop_top_index1)+1:end,:); 
    new_drop_nodal_coordinate2(length(new_drop_top_index2)+1:end,:); 
    new_drop_nodal_coordinate3(length(new_drop_top_index3)+1:end,:); 
    new_drop_nodal_coordinate4(length(new_drop_top_index4)+1:end,:) ;
    drop_nodal_coordinate1(drop_nodal_index,:); drop_nodal_coordinate1(length(drop_top_index1)+1:end,:);
    drop_nodal_coordinate2(drop_nodal_index,:);drop_nodal_coordinate2(length(drop_top_index2)+1:end,:); 
    drop_nodal_coordinate3(drop_nodal_index,:);drop_nodal_coordinate3(length(drop_top_index3)+1:end,:); 
    drop_nodal_coordinate4(drop_nodal_index,:);drop_nodal_coordinate4(length(drop_top_index4)+1:end,:) ;
    col_nodal_coordinate1(length(col_top_index1)+1:end,:); col_nodal_coordinate2(length(col_top_index2)+1:end,:); 
    col_nodal_coordinate3(length(col_top_index3)+1:end,:); col_nodal_coordinate4(length(col_top_index4)+1:end,:) ];
faces = [slab_faces ; new_drop_faces1; new_drop_faces2 ;new_drop_faces3; new_drop_faces4; drop_faces1; drop_faces2 ;drop_faces3; drop_faces4; col_faces1; col_faces2 ;col_faces3; col_faces4];
% faces = [slab_faces ; new_drop_faces1; new_drop_faces2 ;new_drop_faces3; new_drop_faces4; drop_faces1; drop_faces2 ;drop_faces3; drop_faces4;];
nodal_connect = [slab_nodal_connect ;new_drop_nodal_connect1; new_drop_nodal_connect2; new_drop_nodal_connect3; new_drop_nodal_connect4;drop_nodal_connect1; drop_nodal_connect2; drop_nodal_connect3; drop_nodal_connect4;];
%%%
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
function [nodal_connect, faces, bot_node, last_index] = helper4(starting_index,col_x_coord,col_y_coord, col_z_coord, top_index,bot_node)
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

function col_top_index = helper10(drop_bot_coordinate,col_top_coordinate,drop_bot_node)
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
end





function [nodal_connect, faces, bot_node, last_index] = helper5(starting_index, col_x_coord, col_y_coord, col_z_coord, top_index, bot_node)
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
        if(ii == 1 || ii == 3 || ii == 7 || i == 9)
            mapping = [node_one , top_index((row)*(mesh_meta_data(1)+1) + index +  1), top_index((row+1)*(mesh_meta_data(1)+1) + index +  1) , top_index((row+1)*(mesh_meta_data(1)+1) + index ) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
        else
            mapping = [node_one , node_one + 1, top_index((row+1)*(mesh_meta_data(1)+1) + index +  1) , top_index((row+1)*(mesh_meta_data(1)+1) + index ) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
        end
%         mapping = [node_one , node_one + 1, top_index((row+1)*(mesh_meta_data(1)+1) + index +  1) , top_index((row+1)*(mesh_meta_data(1)+1) + index ) , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1 ];
        last_index = max(mapping(:));
        if layer == mesh_meta_data(3)-1
                bot_node =[bot_node node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
        end
            bot_node = unique(bot_node);

        
    else
            node_one = starting_index + (layer-1)*(mesh_meta_data(1)+1)*(mesh_meta_data(2)+1) + row*(mesh_meta_data(1)+1) + index;
            node_two = node_one + (mesh_meta_data(1)+1)*(mesh_meta_data(2)+1);
            mapping = [node_one, node_two, node_one + mesh_meta_data(1) + 2,  node_one + mesh_meta_data(1) + 1 , node_two, node_two + 1,  node_two + mesh_meta_data(1) + 2, node_two + mesh_meta_data(1) + 1];
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
