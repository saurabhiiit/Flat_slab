function [ nodal_connect, nodal_coordinate, faces, slab_mesh_meta_data, reinforcment_element_center, total_node, slab_x_center_coord,slab_y_center_coord, col_bot_node,slab_side_nodes,slab_top_nodes  ] = createMesh(slab_dimension,col_dimension,slab_divisions,col_divisions,col_top_postion, slab_reinforcment_info, col_reinforcment_info )
%%
%CREATEMESH Summary of this function goes here
%   Detailed explanation goes here

% This function divide the slab and column in the mesh and it also consider reinforcment.

% bar positions in the z direction.
col_stirup_spacing = 200;
slab_clear_cover_z = 25;
    

%%% Slab bars positions  in x  y and z direction.
% Total number of bars in y and z direction
% slab_num_bars = floor((slab_dimension(1:2) - 2*slab_reinforcment_info(:,3).')./slab_reinforcment_info(:,2).');

% Total number of bars in y and z direction
% slab_num_bars = floor((slab_dimension(1:2) - 2*slab_reinforcment_info(:,3).')./slab_reinforcment_info(:,2).');

% Update the side covers so as we have equal spacing both the sides.
% temp = (slab_dimension(1:2) -  slab_num_bars.*slab_reinforcment_info(:,2).')/2;
% slab_reinforcment_info(:, 3) = temp;

% slab_bar_poistion_x = (0:slab_num_bars(1))*slab_reinforcment_info(1,2) + slab_reinforcment_info(1,3); % Along x-axis
% slab_bar_poistion_y = (0:slab_num_bars(2))*slab_reinforcment_info(2,2) + slab_reinforcment_info(2,3); % Along y-axis

slab_bar_poistion_x = [slab_reinforcment_info(1,3):slab_reinforcment_info(1,2):slab_dimension(1)/2 ]
slab_bar_poistion_x = sort([slab_bar_poistion_x slab_dimension(1)-slab_reinforcment_info(1,3)-slab_reinforcment_info(1,1):-slab_reinforcment_info(1,2):slab_dimension(1)/2 ])

slab_bar_poistion_y = [slab_reinforcment_info(2,3):slab_reinforcment_info(2,2):slab_dimension(2)/2 ];
slab_bar_poistion_y = sort([slab_bar_poistion_y slab_dimension(2)-slab_reinforcment_info(2,3)-slab_reinforcment_info(2,1):-slab_reinforcment_info(2,2):slab_dimension(2)/2 ]);


% position of bar along the z direction 
if slab_dimension(1)>=slab_dimension(2)
slab_bar_poistion_z = col_top_postion(1,3)+[slab_clear_cover_z slab_clear_cover_z+slab_reinforcment_info(1,1) slab_dimension(3)-slab_clear_cover_z-slab_reinforcment_info(1,1)-slab_reinforcment_info(2,1), slab_dimension(3)-slab_clear_cover_z-slab_reinforcment_info(1,1) ]; % Along z-axis
slab_bar_dia = [slab_reinforcment_info(1,1) slab_reinforcment_info(2,1) slab_reinforcment_info(2,1) slab_reinforcment_info(1,1)];
% slab_tot_bar_postion = unique([slab_bars_poistion_z  ])
else
slab_bars_poistion_z = col_top_postion(1,3)+[slab_clear_cover_z slab_clear_cover_z+slab_reinforcment_info(2,1) slab_dimension(3)-slab_clear_cover_z-slab_reinforcment_info(2,1)-slab_reinforcment_info(1,1), slab_dimension(3)-slab_clear_cover_z-slab_reinforcment_info(2,1) ]; % Along z-axis
slab_bar_dia = [slab_reinforcment_info(2,1) slab_reinforcment_info(1,1) slab_reinforcment_info(1,1) slab_reinforcment_info(2,1)];
end    

% This is done so that to have meshing in x y z direction of the  reinforcment of slab 
temp1 = (slab_bar_poistion_x>(col_top_postion(1,1)-(slab_dimension(1)/12))&slab_bar_poistion_x<(col_top_postion(1,1)-(col_dimension(1)/2)))|(slab_bar_poistion_x>(col_top_postion(1,1)+(col_dimension(1)/2))&slab_bar_poistion_x<(col_top_postion(1,1)+(slab_dimension(1)/12)));
temp2 = (slab_bar_poistion_x>(col_top_postion(2,1)-(slab_dimension(1)/12))&slab_bar_poistion_x<(col_top_postion(2,1)-(col_dimension(1)/2)))|(slab_bar_poistion_x>(col_top_postion(2,1)+(col_dimension(1)/2))&slab_bar_poistion_x<(col_top_postion(2,1)+(slab_dimension(1)/12)));

% temp2 = slab_bar_poistion_x>(col_top_postion(2,1)-(slab_dimension(1)/12))&slab_bar_poistion_x<(col_top_postion(2,1)+(slab_dimension(1)/12));
slab_bar_x = helper6(slab_reinforcment_info(1,1),slab_bar_poistion_x(temp1|temp2)); 
% temp = slab_bar_poistion_x(slab_bar_poistion_x<(col_top_postion(1,1)-(slab_dimension(1)/12))|slab_bar_poistion_x>(col_top_postion(2,1)+(slab_dimension(1)/12))|(slab_bar_poistion_x>(col_top_postion(1,1)+(slab_dimension(1)/12))&slab_bar_poistion_x<(col_top_postion(2,1)-(slab_dimension(1)/12))));
temp = slab_bar_poistion_x(~(temp1|temp2));
slab_bar_x = unique([slab_bar_x temp slab_reinforcment_info(1,1)+temp]) ;


temp1 = (slab_bar_poistion_y>(col_top_postion(1,2)-(slab_dimension(2)/12))&slab_bar_poistion_y<(col_top_postion(1,2)-(col_dimension(2)/2)))|(slab_bar_poistion_y>(col_top_postion(1,2)+(col_dimension(2)/2))&slab_bar_poistion_y<(col_top_postion(1,2)+(slab_dimension(2)/12)));
temp2 = (slab_bar_poistion_y>(col_top_postion(3,2)-(slab_dimension(2)/12))&slab_bar_poistion_y<(col_top_postion(3,2)-(col_dimension(2)/2)))|(slab_bar_poistion_y>(col_top_postion(3,2)+(col_dimension(2)/2))&slab_bar_poistion_y<(col_top_postion(3,2)+(slab_dimension(2)/12)));

% temp1 = slab_bar_poistion_y>(col_top_postion(1,2)-(slab_dimension(2)/12))&slab_bar_poistion_y<(col_top_postion(1,2)+(slab_dimension(2)/12));
% temp2 = slab_bar_poistion_y>(col_top_postion(3,2)-(slab_dimension(2)/12))&slab_bar_poistion_y<(col_top_postion(3,2)+(slab_dimension(2)/12));
slab_bar_y = helper6(slab_reinforcment_info(2,1),slab_bar_poistion_y(temp1|temp2)); 
% temp = slab_bar_poistion_y(slab_bar_poistion_y<(col_top_postion(1,2)-(slab_dimension(2)/12))|slab_bar_poistion_y>(col_top_postion(3,2)+(slab_dimension(2)/12))|(slab_bar_poistion_y>(col_top_postion(1,2)+(slab_dimension(2)/12))&slab_bar_poistion_y<(col_top_postion(3,2)-(slab_dimension(2)/12))));
temp = slab_bar_poistion_y(~(temp1|temp2));
slab_bar_y = unique([slab_bar_y temp slab_reinforcment_info(2,1)+temp]) ;

% slab_bar_y = helper6(slab_reinforcment_info(2,1),slab_bar_poistion_y(slab_bar_poistion_y<=(slab_dimension(2)/3)|slab_bar_poistion_y>=(slab_dimension(2)-(slab_dimension(2)/3))) ); 
% temp = slab_bar_poistion_y(slab_bar_poistion_y>(slab_dimension(2)/3)&slab_bar_poistion_y<(slab_dimension(2)-(slab_dimension(2)/3)));
% slab_bar_y = unique([slab_bar_y temp slab_reinforcment_info(2,1)+temp]) ;

slab_bar_z = helper6(slab_bar_dia,slab_bar_poistion_z ); 

%%% Slab mesh division .
slab_mesh_size = 50 * floor((slab_dimension./slab_divisions)/50); 

slab_x_coord = helper(slab_mesh_size(1), slab_dimension(1),0);
slab_y_coord = helper(slab_mesh_size(2), slab_dimension(2),0);
slab_z_coord = helper(slab_mesh_size(3), slab_dimension(3), col_top_postion(1,3));

%%% Column bars positions.

% For the column 1 
col_stirup_postion_x1 = [(col_top_postion(1,1) - (col_dimension(1)/2))+col_reinforcment_info(2,2) (col_top_postion(1,1)-(col_dimension(1)/2))+col_dimension(1)-col_reinforcment_info(2,2)-col_reinforcment_info(2,1) ];
col_stirup_postion_y1 = [(col_top_postion(1,2) - (col_dimension(2)/2))+col_reinforcment_info(1,2) (col_top_postion(1,2) - (col_dimension(2)/2))+col_dimension(2)-col_reinforcment_info(1,2)-col_reinforcment_info(1,1) ];
% for z direction 
col_num_stirup = floor((col_dimension(3) - 2*col_stirup_spacing.')./col_stirup_spacing.');
col_stirup_postion_z = (0:col_num_stirup)*col_stirup_spacing + col_stirup_spacing; % Along z-axis

col_bar_postion_x1 = [col_stirup_postion_x1(1)+col_reinforcment_info(2,1)  col_stirup_postion_x1(2)-col_reinforcment_info(3,1)];
col_bar_postion_y1 = [col_stirup_postion_y1(1)+col_reinforcment_info(1,1)  col_stirup_postion_y1(2)-col_reinforcment_info(3,1)];
    
% For the column 4 
col_stirup_postion_x4 = [(col_top_postion(4,1) - (col_dimension(1)/2))+col_reinforcment_info(2,2) (col_top_postion(4,1) - (col_dimension(1)/2))+col_dimension(1)-col_reinforcment_info(2,2)-col_reinforcment_info(2,1) ]; 
col_stirup_postion_y4 = [(col_top_postion(4,2) - (col_dimension(2)/2))+col_reinforcment_info(1,2) (col_top_postion(4,2) - (col_dimension(2)/2))+col_dimension(2)-col_reinforcment_info(1,2)-col_reinforcment_info(1,1) ];

col_bar_postion_x4 = [col_stirup_postion_x4(1)+col_reinforcment_info(2,1)  col_stirup_postion_x4(2)-col_reinforcment_info(3,1)];
col_bar_postion_y4 = [col_stirup_postion_y4(1)+col_reinforcment_info(1,1)  col_stirup_postion_y4(2)-col_reinforcment_info(3,1)];

% col total bars divisions
col_tot_bar_postion_x1 = unique([col_stirup_postion_x1 col_reinforcment_info(2,1)+col_stirup_postion_x1 col_bar_postion_x1 col_reinforcment_info(3,1)+col_bar_postion_x1 ]);
col_tot_bar_postion_x2 = unique([col_stirup_postion_x4 col_reinforcment_info(2,1)+col_stirup_postion_x4 col_bar_postion_x4 col_reinforcment_info(3,1)+col_bar_postion_x4 ]);
col_tot_bar_postion_x3 = unique([col_stirup_postion_x1 col_reinforcment_info(2,1)+col_stirup_postion_x1 col_bar_postion_x1 col_reinforcment_info(3,1)+col_bar_postion_x1 ]);
col_tot_bar_postion_x4 = unique([col_stirup_postion_x4 col_reinforcment_info(2,1)+col_stirup_postion_x4 col_bar_postion_x4 col_reinforcment_info(3,1)+col_bar_postion_x4 ]);

col_tot_bar_postion_y1 = unique([col_stirup_postion_y1 col_reinforcment_info(1,1)+col_stirup_postion_y1 col_bar_postion_y1 col_reinforcment_info(3,1)+col_bar_postion_y1 ]);
col_tot_bar_postion_y3 = unique([col_stirup_postion_y4 col_reinforcment_info(1,1)+col_stirup_postion_y4 col_bar_postion_y4 col_reinforcment_info(3,1)+col_bar_postion_y4 ]);
col_tot_bar_postion_y2 = unique([col_stirup_postion_y1 col_reinforcment_info(1,1)+col_stirup_postion_y1 col_bar_postion_y1 col_reinforcment_info(3,1)+col_bar_postion_y1 ]);
col_tot_bar_postion_y4 = unique([col_stirup_postion_y4 col_reinforcment_info(1,1)+col_stirup_postion_y4 col_bar_postion_y4 col_reinforcment_info(3,1)+col_bar_postion_y4 ]);

col_tot_bar_postion_z = unique([col_stirup_postion_z col_reinforcment_info(1,1)+col_stirup_postion_z]);


%%% column mesh division

col_mesh_size = 50 * floor((col_dimension./col_divisions)/50); 

% col 1
col_x_coord1 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(1,1) - (col_dimension(1)/2)));
col_y_coord1 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(1,2) - (col_dimension(2)/2));  
col_z_coord1 = helper(col_mesh_size(3), col_dimension(3), 0);

% col2
col_x_coord2 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(2,1) - (col_dimension(1)/2)));
col_y_coord2 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(2,2) - (col_dimension(2)/2));

col_x_coord3 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(3,1) - (col_dimension(1)/2)));
col_y_coord3 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(3,2) - (col_dimension(2)/2));

col_x_coord4 = helper(col_mesh_size(1), col_dimension(1),(col_top_postion(4,1) - (col_dimension(1)/2)));
col_y_coord4 = helper(col_mesh_size(2), col_dimension(2),col_top_postion(4,2) - (col_dimension(2)/2));

%%% Combining the divisons in the column to column coord to form the combined matrix

col_x_coord1 = unique([col_x_coord1 col_tot_bar_postion_x1 (slab_bar_x(slab_bar_x>(col_top_postion(1,1) - (col_dimension(1)/2))&slab_bar_x<(col_top_postion(1,1) + (col_dimension(1)/2))))]);
col_y_coord1 = unique([col_y_coord1 col_tot_bar_postion_y1 (slab_bar_y(slab_bar_y>(col_top_postion(1,2) - (col_dimension(2)/2))&slab_bar_y<(col_top_postion(1,2) + (col_dimension(2)/2))))]);
col_z_coord1 = unique([col_tot_bar_postion_z, col_z_coord1]);
col_z_coord1 = fliplr(col_z_coord1);

col_nodal_coordinate1 = combvec(col_x_coord1, col_y_coord1, col_z_coord1).';

col_x_coord2 = unique([col_x_coord2 col_tot_bar_postion_x2 (slab_bar_x(slab_bar_x>(col_top_postion(2,1) - (col_dimension(1)/2))&slab_bar_x<(col_top_postion(2,1) + (col_dimension(1)/2))))]);
col_y_coord2 = unique([col_y_coord2 col_tot_bar_postion_y2 (slab_bar_y(slab_bar_y>(col_top_postion(2,2) - (col_dimension(2)/2))&slab_bar_y<(col_top_postion(2,2) + (col_dimension(2)/2))))]);

col_nodal_coordinate2 = combvec(col_x_coord2, col_y_coord2, col_z_coord1).';

col_x_coord3 = unique([col_x_coord3 col_tot_bar_postion_x3 (slab_bar_x(slab_bar_x>(col_top_postion(3,1) - (col_dimension(1)/2))&slab_bar_x<(col_top_postion(3,1) + (col_dimension(1)/2))))]);
col_y_coord3 = unique([col_y_coord3 col_tot_bar_postion_y3 (slab_bar_y(slab_bar_y>(col_top_postion(3,2) - (col_dimension(2)/2))&slab_bar_y<(col_top_postion(3,2) + (col_dimension(2)/2))))]);

col_nodal_coordinate3 = combvec(col_x_coord3, col_y_coord3, col_z_coord1).';

col_x_coord4 = unique([col_x_coord4 col_tot_bar_postion_x4 (slab_bar_x(slab_bar_x>(col_top_postion(4,1) - (col_dimension(1)/2))&slab_bar_x<(col_top_postion(4,1) + (col_dimension(1)/2))))]);
col_y_coord4 = unique([col_y_coord4 col_tot_bar_postion_y4 (slab_bar_y(slab_bar_y>(col_top_postion(4,2) - (col_dimension(2)/2))&slab_bar_y<(col_top_postion(4,2) + (col_dimension(2)/2))))]);

%%%
col_nodal_coordinate4 = combvec(col_x_coord4, col_y_coord4, col_z_coord1).';

%%% Combining the divisons in the slab coord to form the combined matrix

slab_x_coord  = unique([slab_x_coord((slab_x_coord>(col_top_postion(1,1)+(col_dimension(1)/2))&slab_x_coord<(col_top_postion(2,1)-(col_dimension(1)/2)))|slab_x_coord<(col_top_postion(1,1) - (col_dimension(1)/2))|slab_x_coord>(col_top_postion(2,1) + (col_dimension(1)/2))) col_x_coord1 col_x_coord2 slab_bar_x]);
slab_y_coord  = unique([slab_y_coord((slab_y_coord>(col_top_postion(1,2)+(col_dimension(2)/2))&slab_y_coord<(col_top_postion(3,2)-(col_dimension(2)/2)))|slab_y_coord<(col_top_postion(1,1) - (col_dimension(2)/2))|slab_y_coord>(col_top_postion(3,2) + (col_dimension(2)/2))) col_y_coord1 col_y_coord3 slab_bar_y]);
slab_z_coord = unique([slab_bar_z, slab_z_coord]);

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
  
%%% Top index of the column
col_top_index1 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate1);
col_top_index2 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate2);
col_top_index3 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate3);
col_top_index4 = helper3(slab_x_coord, slab_y_coord, slab_nodal_coordinate, col_top_coordinate4);

%%% Nodal connectivity of the each four column
col_bot_node=[];
[col_nodal_connect1, col_faces1, col_bot_node1, col_last_index1] = helper4( max(slab_nodal_connect(:)), col_x_coord1, col_y_coord1, col_z_coord1, col_top_index1,col_bot_node);
[col_nodal_connect2, col_faces2, col_bot_node2, col_last_index2] = helper4( col_last_index1, col_x_coord2, col_y_coord2, col_z_coord1, col_top_index2,col_bot_node);
[col_nodal_connect3, col_faces3, col_bot_node3, col_last_index3] = helper4( col_last_index2, col_x_coord3, col_y_coord3, col_z_coord1, col_top_index3,col_bot_node);
[col_nodal_connect4, col_faces4, col_bot_node4, col_last_index4] = helper4( col_last_index3, col_x_coord4, col_y_coord4, col_z_coord1, col_top_index4,col_bot_node);
col_bot_node = [col_bot_node1 col_bot_node2 col_bot_node3 col_bot_node4];
% disp(col_bot_node3);
bot_node_len = [length(col_bot_node1) length(col_bot_node2) length(col_bot_node3) length(col_bot_node4)]

%%%

%%% Combining the nodal coordinate of slab and all columns to one
nodal_coordinate = [slab_nodal_coordinate ;col_nodal_coordinate1(length(col_top_index1)+1:end,:); col_nodal_coordinate2(length(col_top_index2)+1:end,:); col_nodal_coordinate3(length(col_top_index3)+1:end,:); col_nodal_coordinate4(length(col_top_index4)+1:end,:) ];
faces = [slab_faces ; col_faces1; col_faces2 ;col_faces3; col_faces4];
nodal_connect = [slab_nodal_connect ;col_nodal_connect1; col_nodal_connect2; col_nodal_connect3; col_nodal_connect4];
%%%

% length(slab_nodal_connect)
% length(nodal_connect)
%%% Calculating the middle point of the element having reinforcment

% Slab Reinforcment marking
slab_reinforcement_center_x =  slab_bar_poistion_x + (slab_reinforcment_info(1,1))/2;
[slab_center_x_x] = helper8( slab_reinforcement_center_x, slab_reinforcment_info(1,1), slab_x_coord );
slab_center_x_y = (slab_y_coord(3:end-1) + slab_y_coord(2:end-2))/2;

slab_reinforcement_center_y =  slab_bar_poistion_y + (slab_reinforcment_info(2,1))/2;
[slab_center_y_y] = helper8( slab_reinforcement_center_y, slab_reinforcment_info(2,1), slab_y_coord );
slab_center_y_x = (slab_x_coord(3:end-1) + slab_x_coord(2:end-2))/2;

if slab_dimension(1)>=slab_dimension(2)
    slab_center_x_z = [(slab_z_coord(2)+slab_z_coord(3))/2  (slab_z_coord(end-1)+slab_z_coord(end-2))/2];
    slab_center_y_z = [(slab_z_coord(3)+slab_z_coord(4))/2  (slab_z_coord(end-2)+slab_z_coord(end-3))/2];
else
    slab_center_y_z = [(slab_z_coord(2)+slab_z_coord(3))/2  (slab_z_coord(end-1)+slab_z_coord(end-2))/2] ;
    slab_center_x_z = [(slab_z_coord(3)+slab_z_coord(4))/2  (slab_z_coord(end-2)+slab_z_coord(end-3))/2];
end

slab_x_center_coord = [combvec(unique(slab_center_x_x), unique(slab_center_x_y),unique(slab_center_x_z(1)))' ; combvec(slab_center_x_x(slab_center_x_x<slab_dimension(1)/3|slab_center_x_x>(slab_dimension(1)-(slab_dimension(1)/3))), slab_center_x_y,slab_center_x_z(2))'];
slab_y_center_coord = [combvec(slab_center_y_x, slab_center_y_y,slab_center_y_z(1))' ; combvec(slab_center_y_x,slab_center_y_y(slab_center_y_y<slab_dimension(2)/3|slab_center_y_y>(slab_dimension(2)-(slab_dimension(2)/3))),slab_center_y_z(2))' ];

% disp(unique(slab_center_x_x))
% conserdering the slabs elements in the sequence
if slab_dimension(1)>=slab_dimension(1)
    slab_elements = [combvec(unique(slab_center_x_x), unique(slab_center_x_y),unique(slab_center_x_z(1)))' ;combvec(unique(slab_center_y_x), unique(slab_center_y_y),unique(slab_center_y_z(1)))' ; 
    combvec(unique(slab_center_y_x),unique(slab_center_y_y(slab_center_y_y<slab_dimension(2)/3|slab_center_y_y>(slab_dimension(2)-(slab_dimension(2)/3)))),unique(slab_center_y_z(2)))' ;combvec(unique(slab_center_x_x(slab_center_x_x<slab_dimension(1)/3|slab_center_x_x>(slab_dimension(1)-(slab_dimension(1)/3)))), unique(slab_center_x_y),unique(slab_center_x_z(2)))'];
else
    slab_elements=[combvec(unique(slab_center_y_x), unique(slab_center_y_y),unique(slab_center_y_z(1)))' ; combvec(unique(slab_center_x_x), unique(slab_center_x_y),unique(slab_center_x_z(1)))' ; 
    combvec(unique(slab_center_x_x(slab_center_x_x<slab_dimension(1)/3|slab_center_x_x>(slab_dimension(1)-(slab_dimension(1)/3)))), unique(slab_center_x_y,slab_center_x_z(2)))'; combvec(unique(slab_center_y_x),unique(slab_center_y_y(slab_center_y_y<slab_dimension(2)/3|slab_center_y_y>(slab_dimension(2)-(slab_dimension(2)/3)))),unique(slab_center_y_z(2)))' ];
end

% col Reinforcment marking
% bar elements
% disp(col_bar_postion_x1 )
col_reinforcement_center_bar_x1 = col_bar_postion_x1 + (col_reinforcment_info(3,1))/2;
[col_bar_center_x1] = helper8( col_reinforcement_center_bar_x1, col_reinforcment_info(3,1), col_x_coord1 );
 
col_reinforcement_center_bar_y1 = col_bar_postion_y1 + col_reinforcment_info(3,1)/2;
[col_bar_center_y1] = helper8( col_reinforcement_center_bar_y1, col_reinforcment_info(3,1), col_y_coord1 );

col_reinforcement_center_bar_x4 = col_bar_postion_x4 + col_reinforcment_info(3,1)/2;
[col_bar_center_x4] = helper8( col_reinforcement_center_bar_x4, col_reinforcment_info(3,1), col_x_coord4 );

col_reinforcement_center_bar_y4 = col_bar_postion_y4 + col_reinforcment_info(3,1)/2;
[col_bar_center_y4] = helper8( col_reinforcement_center_bar_y4, col_reinforcment_info(3,1), col_y_coord4 );

col_bar_center_z = (col_z_coord1(1:end-1) + col_z_coord1(2:end))/2;



% stirup elements 
%col1
%x perpend bar
col_reinforcement_center_stirup_x1 = col_stirup_postion_x1 + (col_reinforcment_info(2,1)/2);
[col_stirup_center_x_x1] = helper8( col_reinforcement_center_stirup_x1, col_reinforcment_info(2,1), col_x_coord1 );
% col_y_coord1
col_stirup_center_x_y1_temp = col_y_coord1(col_y_coord1>=col_top_postion(1,2)-(col_dimension(2)/2)+col_reinforcment_info(1,2)&col_y_coord1<=col_top_postion(1,2)+(col_dimension(2)/2)-col_reinforcment_info(1,2));
col_stirup_center_x_y1 = (col_stirup_center_x_y1_temp(1:end-1) + col_stirup_center_x_y1_temp(2:end))/2;

%y perpend bar
col_reinforcement_center_stirup_y1 = col_stirup_postion_y1 + (col_reinforcment_info(1,1)/2);
[col_stirup_center_y_y1] = helper8( col_reinforcement_center_stirup_y1, col_reinforcment_info(1,1), col_y_coord1 );
col_stirup_center_y_x1_temp = col_x_coord1(col_x_coord1>=col_top_postion(1,1)-(col_dimension(1)/2)+col_reinforcment_info(2,2)&col_x_coord1<=col_top_postion(1,1)+(col_dimension(1)/2)-col_reinforcment_info(2,2));
col_stirup_center_y_x1 = (col_stirup_center_y_x1_temp(1:end-1) + col_stirup_center_y_x1_temp(2:end))/2;


%col4
%x perpend bar
col_reinforcement_center_stirup_x4 = col_stirup_postion_x4 + (col_reinforcment_info(2,1)/2);
[col_stirup_center_x_x4] = helper8( col_reinforcement_center_stirup_x4, col_reinforcment_info(2,1), col_x_coord4 );
col_stirup_center_x_y4_temp = col_y_coord4(col_y_coord4>=col_top_postion(4,2)-(col_dimension(2)/2)+col_reinforcment_info(1,2)&col_y_coord4<=col_top_postion(4,2)+(col_dimension(1)/2)-col_reinforcment_info(1,2));
col_stirup_center_x_y4 = (col_stirup_center_x_y4_temp(1:end-1) + col_stirup_center_x_y4_temp(2:end))/2;

%y perpend bar
col_reinforcement_center_stirup_y4 = col_stirup_postion_y4 + (col_reinforcment_info(1,1)/2);
[col_stirup_center_y_y4] = helper8( col_reinforcement_center_stirup_y4, col_reinforcment_info(1,1), col_y_coord4 );
col_stirup_center_y_x4_temp = col_x_coord4(col_x_coord4>=col_top_postion(4,1)-(col_dimension(1)/2)+col_reinforcment_info(2,2)&col_x_coord4<=col_top_postion(4,1)+(col_dimension(1)/2)-col_reinforcment_info(2,2));
col_stirup_center_y_x4 = (col_stirup_center_y_x4_temp(1:end-1) + col_stirup_center_y_x4_temp(2:end))/2;

col_stirup_center_z = col_stirup_postion_z + col_reinforcment_info(1,1)/2;

%col1
% unique(col_stirup_center_x_x1)


col_elements1  = [combvec(unique(col_bar_center_x1),unique(col_bar_center_y1),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x1),unique(col_stirup_center_x_y1),fliplr(unique(col_stirup_center_z)))';
combvec(unique(col_stirup_center_y_x1),unique(col_stirup_center_y_y1),fliplr(unique(col_stirup_center_z)))'];
    
col_elements2  = [combvec(unique(col_bar_center_x4),unique(col_bar_center_y1),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x4),unique(col_stirup_center_x_y1),fliplr(unique(col_stirup_center_z)))';
combvec(unique(col_stirup_center_y_x4),unique(col_stirup_center_y_y1),fliplr(unique(col_stirup_center_z)))'];

col_elements3  = [combvec(unique(col_bar_center_x1),unique(col_bar_center_y4),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x1),unique(col_stirup_center_x_y4),fliplr(unique(col_stirup_center_z)))';
combvec(unique(col_stirup_center_y_x1),unique(col_stirup_center_y_y4),fliplr(unique(col_stirup_center_z)))'];

col_elements4  = [combvec(unique(col_bar_center_x4),unique(col_bar_center_y4),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x4),unique(col_stirup_center_x_y4),fliplr(unique(col_stirup_center_z)))';
combvec(unique(col_stirup_center_y_x4),unique(col_stirup_center_y_y4),fliplr(unique(col_stirup_center_z)))'];


% col_elements1  = [combvec(unique(col_reinforcement_center_bar_x1),unique(col_reinforcement_center_bar_y1),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x1),unique(col_stirup_center_x_y1),fliplr(unique(col_stirup_center_z)))';
% combvec(unique(col_stirup_center_y_x1),unique(col_stirup_center_y_y1),fliplr(unique(col_stirup_center_z)))'];
    
% col_elements2  = [combvec(unique(col_reinforcement_center_bar_x4),unique(col_reinforcement_center_bar_y1),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x4),unique(col_stirup_center_x_y1),fliplr(unique(col_stirup_center_z)))';
% combvec(unique(col_stirup_center_y_x4),unique(col_stirup_center_y_y1),fliplr(unique(col_stirup_center_z)))'];

% col_elements3  = [combvec(unique(col_reinforcement_center_bar_x1),unique(col_reinforcement_center_bar_y4),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x1),unique(col_stirup_center_x_y4),fliplr(unique(col_stirup_center_z)))';
% combvec(unique(col_stirup_center_y_x1),unique(col_stirup_center_y_y4),fliplr(unique(col_stirup_center_z)))'];

% col_elements4  = [combvec(unique(col_reinforcement_center_bar_x4),unique(col_reinforcement_center_bar_y4),unique(col_bar_center_z))';combvec(unique(col_stirup_center_x_x4),unique(col_stirup_center_x_y4),fliplr(unique(col_stirup_center_z)))';
% combvec(unique(col_stirup_center_y_x4),unique(col_stirup_center_y_y4),fliplr(unique(col_stirup_center_z)))'];


reinforcment_element_center = [slab_elements;col_elements1;col_elements2;col_elements3;col_elements4];
total_node = [length(slab_nodal_connect) length(slab_elements);length(col_nodal_connect1) length(col_elements1); length(col_nodal_connect2) length(col_elements2); length(col_nodal_connect3) length(col_elements3); length(col_nodal_connect4) length(col_elements4)];

end



function [coordinates] = helper(mesh_size, dimension,coordinates)
    flag = 0;
    while(flag == 0)
        coord_val = min(coordinates(end) + mesh_size , coordinates(1) + dimension);
        coordinates(end + 1) = coord_val;
        flag = coordinates(end) >= (coordinates(1) + dimension);
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



function [coordinates] = helper6(bar_dia, bar_pos)
    coordinates = [];
    temp = -1:1:2;
    if length(bar_dia)==1
        for i=1:length(bar_pos)
            coordinates =  [coordinates, (bar_pos(i) + bar_dia*temp)];
        end
    else  
        for i=1:length(bar_pos)
            coordinates =  [coordinates, bar_pos(i),(bar_pos(i) + bar_dia(i))];
        end
        coordinates = unique(coordinates);
    end
end

function [slab_nodal_connect , slab_faces, slab_mesh_meta_data,slab_side_nodes,slab_top_nodes] = helper7(slab_x_coord, slab_y_coord, slab_z_coord)
slab_no_elements = (length(slab_x_coord) - 1) * (length(slab_y_coord) - 1) * (length(slab_z_coord) - 1);
disp(slab_no_elements)

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

function [slab_center_x_final] = helper8( slab_reinforcement_center_x, slab_reinforcment_info, slab_x_coord )
slab_center_x_final = [];
for ii = 1:length(slab_reinforcement_center_x)
    if (length(slab_x_coord(slab_x_coord>slab_reinforcement_center_x(ii)-slab_reinforcment_info/2&slab_x_coord<slab_reinforcement_center_x(ii)-slab_reinforcment_info/2)) >0)
        slab_center_x = slab_x_coord(slab_x_coord>=slab_reinforcement_center_x(ii)-slab_reinforcment_info/2&slab_x_coord<=slab_reinforcement_center_x(ii)-slab_reinforcment_info/2);
        slab_center_x_final = [slab_center_x_final slab_center_x(1:end-1)+((slab_center_x(2:end)-slab_center_x(1:end-1))/2)];    
    else
        slab_center_x_final = [slab_center_x_final slab_reinforcement_center_x];
    end
end
end

