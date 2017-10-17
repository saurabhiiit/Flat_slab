function [ global_stiff ] = stiff_global(global_stiff,nodal_coordinate, nodal_connect,  element_mapping, mod_of_elas )
%STIFF_GLOBAL Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(nodal_connect)

[stiff] = octa_element_stiff(nodal_connect(i,:), nodal_coordinate,mod_of_elas(i));

for j = 1:24
  for k = 1:24
  global_stiff(element_mapping(i,j),element_mapping(i,k)) = global_stiff(element_mapping(i,j),element_mapping(i,k))  + stiff(j,k);
  % disp([num2str(i) , '\t', num2str(j) , '\t', num2str(k)])
  end
end
end
end