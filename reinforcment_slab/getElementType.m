function element_type_steel = getElementType(nodal_coordinate, nodal_connect, reinforcment_element_center, total_node)
%**************************************************************************
% Assigns element type to each elements. Returns a boolean vector storing
% whether the element is steel or not. 0 represents concrete element.
%**************************************************************************
%
% Input parameters:
% nodal_connect      - Nodal connectivity matrix in which each row has the
%                      nodes present in that element.
% nodal_coordinate   - Nodal coordinate matrix having each column as the
%                      x,y and z coordinate value of the repective node. 
% bar_position       - Reinforcement bars position in each directions.
% thickness          - Thickness of the wall used to determine how many
%                      layers of reinforcement is present and their 
%                      orientation.
% Output:
% element_type_steel - A boolean vector storing whether the element is 
%                      steel or not. 0 represents concrete element.

no_elements = length(nodal_connect);
element_type_steel = zeros(no_elements, 1);

kk=[];
for ii=1:no_elements
    mid_x = (nodal_coordinate(nodal_connect(ii,1),1)+nodal_coordinate(nodal_connect(ii,2),1))/2;
    mid_y = (nodal_coordinate(nodal_connect(ii,1),2)+nodal_coordinate(nodal_connect(ii,3),2))/2;
    mid_z = (nodal_coordinate(nodal_connect(ii,1),3)+nodal_coordinate(nodal_connect(ii,5),3))/2;

    if ~sum(reinforcment_element_center(:,1)==mid_x&reinforcment_element_center(:,2)==mid_y&reinforcment_element_center(:,3)==mid_z)==0
         % disp(ii)
        element_type_steel(ii)=1;
        % kk = [kk ii];
    % if sum(reinforcment_element_center(:,1)==mid_x&reinforcment_element_center(:,2)==mid_y&reinforcment_element_center(:,3)==mid_z)>1
    %     find(reinforcment_element_center(:,1)==mid_x&reinforcment_element_center(:,2)==mid_y&reinforcment_element_center(:,3)==mid_z)
    % end
    else 
    end
end        
% disp(length(kk));
% disp(length(nonzeros(element_type_steel)));
end

% k=1;
% size_temp = total_node(1,2)+total_node(2,2)+total_node(3,2)+total_node(4,2)+total_node(5,2);
% index_steel = zeros(size_temp,1);
% for ii = 1:no_elements
% % for ii = 1:1
%     first_point_x =  nodal_coordinate(nodal_connect(ii, 1), 1);
%     second_point_x =  nodal_coordinate(nodal_connect(ii, 2), 1);
%     first_point_y =  nodal_coordinate(nodal_connect(ii, 1), 2);
%     second_point_y =  nodal_coordinate(nodal_connect(ii, 4), 2);
%     first_point_z =  nodal_coordinate(nodal_connect(ii, 1), 3);
%     second_point_z =  nodal_coordinate(nodal_connect(ii, 5), 3);
%     for k=1:total_node(1,2)+total_node(2,2)+total_node(3,2)+total_node(4,2)+total_node(5,2)
%     if(first_point_x<reinforcment_element_center(k,1)&reinforcment_element_center(k,1)<second_point_x& first_point_y<reinforcment_element_center(k,2)&reinforcment_element_center(k,2)<second_point_y&first_point_z<reinforcment_element_center(k,3)&reinforcment_element_center(k)<second_point_z)
%         element_type_steel(ii) = 1;
%         % k = k +1;
%          % disp(ii)
%          index_steel(k)=1;
%          break;
%     else
%         element_type_steel(ii) = 0;
%     end
% end

% end
% end





% k=1;
% for ii = 1:total_node(1,1)
% % for ii = 1:1
%     first_point_x =  nodal_coordinate(nodal_connect(ii, 1), 1);
%     second_point_x =  nodal_coordinate(nodal_connect(ii, 2), 1);
%     first_point_y =  nodal_coordinate(nodal_connect(ii, 1), 2);
%     second_point_y =  nodal_coordinate(nodal_connect(ii, 4), 2);
%     first_point_z =  nodal_coordinate(nodal_connect(ii, 1), 3);
%     second_point_z =  nodal_coordinate(nodal_connect(ii, 5), 3);
%     if(first_point_x<reinforcment_element_center(k,1)&reinforcment_element_center(k,1)<second_point_x& first_point_y<reinforcment_element_center(k,2)&reinforcment_element_center(k,2)<second_point_y&first_point_z<reinforcment_element_center(k,3)&reinforcment_element_center(k)<second_point_z)
%         element_type_steel(ii) = 1;
%         k = k +1;
%          disp(ii)
%     else
%         element_type_steel(ii) = 0;
%     end
% end

% for ij = 1:4
% element_type_steel = helper9(element_type_steel,nodal_connect,nodal_coordinate,total_node(ij,:),total_node(ij+1,:),reinforcment_element_center);
% end
% end

% function [element_type_steel] = helper9(element_type_steel,nodal_connect,nodal_coordinate,total_node_start,total_node_end,reinforcment_element_center)
% for ii1 = total_node_start(1)+1:total_node_end(1)

%     first_point_x =  nodal_coordinate(nodal_connect(ii1, 1), :);
%     second_point_x =  nodal_coordinate(nodal_connect(ii1, 2), :);
%     first_point_y =  nodal_coordinate(nodal_connect(ii1, 1), :);
%     second_point_y =  nodal_coordinate(nodal_connect(ii1, 4), :);
%     first_point_z =  nodal_coordinate(nodal_connect(ii1, 1), :);
%     second_point_z =  nodal_coordinate(nodal_connect(ii1, 5), :);
%     for k = total_node_start(2)+1:total_node_end(2)
%         if(first_point_x<reinforcment_element_center(k,1)&reinforcment_element_center(k,1)<second_point_x& first_point_y<reinforcment_element_center(k,2)&reinforcment_element_center(k,2)<second_point_y&first_point_z<reinforcment_element_center(k,3)&reinforcment_element_center(k)<swcond_point_z)
%             element_type_steel(ii1) = 1;
%             break;
%         else
%             element_type_steel(ii1) = 0;
    
%         end
%     end
% end
% end