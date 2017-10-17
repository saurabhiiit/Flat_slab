%**********Analysis of flat slab using FEM Nonlinear analysis without drop and without reinforcement*******************
clear variables;
clear global;
clc;
   
%**********************Inputs******************************
 
%% Input
choice = 0;
while(~(choice == 1 || choice == 2))
    choice = input('[1]. Provide new input values, or\n[2]. Use default values\nChoose [1/2]:    ');
end
%*******Input of flat slab *******************************
if choice == 1
slabx = input('Enter the length of slab in mm: '); %lx = slabl(x) - (slabb(y)/3)
slaby = input('Enter the breadth of slab in mm: ');%ly = 2*slabb(y)/3
slabz = input('Enter the height of slab in mm: ');

slabsteelE = input('Enter the Youngs modulus of steel in N/mm2: ');
slabpoisratio = input('Enter the poisons ratio: ');

slabxspacing = input('Enter the x direction horizontal spacing in mm: ');
slabxdia = input('Enter the x direction main bar diameter in mm: ');
slabyspacing = input('Enter the y direction horizontal spacing in mm: ');
slabydia = input('Enter the y direction main bar diameter in mm: ');

slabfck = input('Enter the grade of slab concrete: ');
slabfy = input('Enter the grade of slab steel: ');

%*********Input of column************************
noc = input('Enter the number of column: ');
colx = input('Enter the length of the column in mm: ');
coly = input('Enter the breadth of the column in mm: ');
colz = input('Enter the height of the column in mm: ');

colhead = input('Enter y for column head ');
if colhead == 'y'
    colhead = 1;
else
    colhead = 0;
end
droppanel = input('Enter y for drop panel ');
if droppanel == 'y'
    droppanel = 1;
else
    droppanel = 0;
end

colsteelE = input('Enter the Young''s modulus of steel in N/mm2: ');
colpoisratio = input('Enter the poisons ratio: ');

colvspacing = input('Enter the vertical direction main bat spacing in mm: ');
colvdia = input('Enter the vertical main bar diameter in mm: ');
coltspacing = input('Enter the tie bar spacing in mm: ');
coltdia = input('Enter the tie bar diameter in mm: ');

colfck = input('Enter the grade of slab concrete: ');
colfy = input('Enter the grade of slab steel: ');

else
slabx = 7500;
slaby = 7500;
slabz = 200;

slabsteelE = 2 * 10^5;
slabpoisratio = 0.3;
steel_E = 2 * 10^5;
slabxspacing =300;
slabxdia = 6;
slabyspacing = 300;
slabydia = 6;

slabfck = 25;
slabfy = 415;	

%*********Input of column************************
noc = 4;
colx = 500;
coly = 500;
colz = 3000;

colhead = 0;
droppanel = 0;

colsteelE = 2 * 10^5;
colpoisratio = 0.3;

colfck = 25;
colfy = 415;

dropz = 100;
end

conc_yield_strain = 0.002;
conc_E = 5000 * sqrt(slabfck);
% conc_E = 5000 * sqrt(conc_grade);
steel_grade = 415;

steel_Et = steel_E / 5;
conc_Et = conc_E / 10;

steel_yield_strain = steel_grade / steel_E;

%%% Density2
conc_d = 2500/(10^9);
steel_d = 7850/(10^9);
%%%

slab_dimension = [slabx, slaby ,slabz];
drop_dimension = [slabx/6, slaby/6, dropz];
drop_top_postion = [slaby/6 slaby/6 colz+dropz;
                   (slabx - slaby/6) slaby/6 colz+dropz;
                   slaby/6 5*slaby/6 colz+dropz;
                   (slabx-slaby/6) 5*slaby/6 colz+dropz];


col_dimension = [colx, coly, colz];
% col_top_postion = [slaby/6 slaby/6 colz;
% 				   (slabx - slaby/6) slaby/6 colz;
% 				   slaby/6 5*slaby/6 colz;
% 				   (slabx-slaby/6) 5*slaby/6 colz];
col_top_postion = [1250 1250 colz;
				   (slabx - 1250) 1250 colz;
				   1250 slaby-1250 colz;
				   (slabx - 1250) (slaby - 1250) colz];

slab_div_x = 15;
slab_div_y = 15;
slab_div_z = 1; % This should vary according to the bars positioning.
slab_divisions = [slab_div_x, slab_div_y, slab_div_z];

col_div_x = 1;
col_div_y = 1;
col_div_z = 15; % This should vary according to the bars positioning.
col_divisions = [col_div_x, col_div_y, col_div_z];

%Reinforcement

%slab
slab_dia_x = 12; % bar perpendicular to x
slab_dia_y = 12;
slab_side_cover_x = 25;% Veticle bars side cover.
slab_side_cover_y = 25;
slab_spacing_x = 300;
slab_spacing_y = 300;
slab_reinforcment_info = [slab_dia_x , slab_spacing_x, slab_side_cover_x;
                     slab_dia_y , slab_spacing_y, slab_side_cover_y];

%column

col_dia_x = 8;
col_dia_y = 8;
col_dia_z = 25;
col_side_cover_x = 40;
col_side_cover_y = 40;
col_side_cover_z = 0;

col_reinforcment_info = [col_dia_x, col_side_cover_x ; col_dia_y col_side_cover_y ; col_dia_z col_side_cover_z ];
%% MESH 
tic;
disp('Creating Mesh...');
[nodal_connect, nodal_coordinate, faces, slab_mesh_meta_data, reinforcment_element_center, total_node, slab_x_center_coord,slab_y_center_coord,col_bot_node,slab_side_nodes,slab_top_nodes] = createMesh(slab_dimension,col_dimension,slab_divisions,col_divisions,col_top_postion,  slab_reinforcment_info, col_reinforcment_info );
disp('Mesh created');
toc;

tic;
%% Draw MESH
disp('Drawing Mesh...');
 draw3DMesh(nodal_coordinate, faces);
disp('Mesh Drawn');
toc;
tic;
%% Get element type
disp('Getting each element type...');
element_type_steel = getElementType(nodal_coordinate, nodal_connect, reinforcment_element_center, total_node);
disp('Done!');
toc;

%%
% Total number of elements will be equal to the the size of the nodal_coordinate matrix.
no_elements = length(nodal_connect);
% disp(no_elements)

total_no_nodes = length(nodal_coordinate);
% disp(total_no_nodes)

total_nodal_displ = zeros(total_no_nodes*3, 1);
[element_mapping, col_bot_dof] = ElementMapping(nodal_connect, no_elements, col_bot_node);

% Initialize the E vector for each element with elastic modulus of
% elasticity.
element_mod_of_elas = steel_E.*element_type_steel(1:no_elements) + conc_E.*(~element_type_steel(1:no_elements));

%% Stiffness Matrix Calculation
global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);


%%% Applying loads
% Dead Load
% It is acting only on bot surfaces of the each elements.
load = zeros(total_no_nodes*3, 1);
mass_slab = 0;
element_density = steel_d.*element_type_steel(1:no_elements) + conc_d.*(~element_type_steel(1:no_elements));
disp('Calculating the gravity Load...');
tic;
for ii=1:no_elements
    dim_x = (nodal_coordinate(nodal_connect(ii,2),1)-nodal_coordinate(nodal_connect(ii,1),1));
    dim_y = (nodal_coordinate(nodal_connect(ii,3),2)-nodal_coordinate(nodal_connect(ii,1),2));
        dim_z = abs(nodal_coordinate(nodal_connect(ii,5),3)-nodal_coordinate(nodal_connect(ii,1),3));
    load(3*nodal_connect(ii,1)) = load(3*nodal_connect(ii,1))-(element_density(ii)*9.81*dim_x*dim_y*dim_z)/8;      
    load(3*nodal_connect(ii,2)) = load(3*nodal_connect(ii,2))-(element_density(ii)*9.81*dim_x*dim_y*dim_z)/8;      
    load(3*nodal_connect(ii,3)) = load(3*nodal_connect(ii,3))-(element_density(ii)*9.81*dim_x*dim_y*dim_z)/8;      
    load(3*nodal_connect(ii,4)) = load(3*nodal_connect(ii,4))-(element_density(ii)*9.81*dim_x*dim_y*dim_z)/8;
    if ii <= (slab_mesh_meta_data(1))*(slab_mesh_meta_data(2))*(slab_mesh_meta_data(3))
        mass_slab =mass_slab+ (element_density(ii)*dim_x*dim_y*dim_z);      
    end
end

% Live Load
% It is acting only on top surfaces of the top elements of slab
% live_load = 4000/(1*10^6);
% for ii = slab_top_nodes
%     dim_x = (nodal_coordinate(nodal_connect(ii,6),1)-nodal_coordinate(nodal_connect(ii,5),1));
%     dim_y = (nodal_coordinate(nodal_connect(ii,8),2)-nodal_coordinate(nodal_connect(ii,5),2));
%     load(3*nodal_connect(ii,5)) = load(3*nodal_connect(ii,5)) - (live_load*dim_x*dim_y)/4;
%     load(3*nodal_connect(ii,6)) = load(3*nodal_connect(ii,6)) - (live_load*dim_x*dim_y)/4;
%     load(3*nodal_connect(ii,7)) = load(3*nodal_connect(ii,7)) - (live_load*dim_x*dim_y)/4;
%     load(3*nodal_connect(ii,8)) = load(3*nodal_connect(ii,8)) - (live_load*dim_x*dim_y)/4;
% end

% Lateral Load
% 1mg acting along the positive x direction on left x-z plane.
lateral_load_mod  = 0.01;
lateral_load = zeros(total_no_nodes*3, 1);
for ii = slab_side_nodes
    dim_y = (nodal_coordinate(nodal_connect(ii,4),2)-nodal_coordinate(nodal_connect(ii,1),2));
    dim_z = (nodal_coordinate(nodal_connect(ii,5),3)-nodal_coordinate(nodal_connect(ii,1),3));
    lateral_load(3*nodal_connect(ii,1)-2) = lateral_load(3*nodal_connect(ii,1)-2) + (lateral_load_mod*dim_y*dim_z)/4;
    lateral_load(3*nodal_connect(ii,4)-2) = lateral_load(3*nodal_connect(ii,4)-2) + (lateral_load_mod*dim_y*dim_z)/4;
    lateral_load(3*nodal_connect(ii,5)-2) = lateral_load(3*nodal_connect(ii,5)-2) + (lateral_load_mod*dim_y*dim_z)/4;
    lateral_load(3*nodal_connect(ii,8)-2) = lateral_load(3*nodal_connect(ii,8)-2) + (lateral_load_mod*dim_y*dim_z)/4;
end

max_displ = [];
total_max_strain = zeros(1, no_elements);
each_ele_strain = zeros(8, 6, no_elements);
total_max_stress = zeros(1, no_elements);
each_ele_stress = zeros(8, 6, no_elements);
count = 0;
% load = zeros(total_no_nodes*3, 1);
% return;
i=0;
while(count==0)
    i=i+1;
    disp(i);
    if i>1
    load = lateral_load;
    end
    residual_force = load;
    while(residual_force(abs(residual_force) > 0))
        load = residual_force;
        
        %% Stiffness Matrix Calculation
         global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);
        %% Boundary conditions
        disp('Applying boundary conditions...');
        tic
        [global_stiff_bc, load_bc] = boundary_conditions( global_stiff, load, col_bot_dof);
        toc;
        disp('Done!');

        %% Solving linear equation
        disp('Solving for nodal displacement...');
        tic
        % Dont do inverse of global matrix as the resultant matrix may not be a
        % sparse matrix and we can't store such huge dense matrix.
        nodal_displ = global_stiff_bc\load_bc;
        total_nodal_displ = total_nodal_displ + nodal_displ;
%         nodal_displ = total_nodal_displ;
        toc
        disp('Done!');

        %% Finding out stress and strain values for each elements
        disp('Calculating element stresses and strains...');
        tic
        count = 0;
          for ii = 1:no_elements
            ele_nodal_disp = nodal_displ(element_mapping(ii, :));
            [max_ele_strain, ele_strain] = ElementStressStrain(nodal_coordinate(nodal_connect(ii, :).', :), ele_nodal_disp);
%             max_strain = max(max(abs(each_ele_strain(:, :, ii) +ele_strain)));
            total_max_strain(ii) = total_max_strain(ii) + max_ele_strain;
%             total_max_strain(ii) = max_strain;
            % Check if element is going into non lienar state or not. If it is then
            % update the element modulus of elasticity or the stress-strain slope.
            % Also do calcualtions here so that we can find out the force value
            % using stress. This will be used to find out the residual force.
            if(total_max_strain(ii) > conc_yield_strain )
                element_mod_of_elas(ii) = conc_Et;
%                 ii
                count = count + 1;
            end
            each_ele_stress(:, :, ii) = each_ele_stress(:, :, ii) + (cal_D(element_mod_of_elas(ii))*ele_strain')';
            total_max_stress(ii) = total_max_stress(ii)+ max(max(abs(cal_D(element_mod_of_elas(ii))*ele_strain')));

            each_ele_strain(:, :, ii) = each_ele_strain(:, :, ii) + ele_strain;
        end      
        toc
        disp('Done!');
        fprintf('Number of elements that went into non-linear state in this iteration: %d\n', count);
        
        %% Calculating Residual Force
        global_stiff = getGlobalStiff(nodal_coordinate, nodal_connect, element_mod_of_elas);
        disp('Applying boundary conditions...');
        tic
        [global_stiff_, load_] = boundary_conditions(global_stiff, load, col_bot_dof);
        toc;
        disp('Done!');
        internal_force = global_stiff_*nodal_displ;
        residual_force = load_bc - internal_force;
        % Remove any floating point error load values
        residual_force(abs(residual_force)<1e-6) = 0;
        % 
%         max(total_nodal_displ)
        
        
    end
    max_displ(end+1) = max(total_nodal_displ);
end
nodal_displ = total_nodal_displ;

%% Desplaying results

nodal_delta_x = nodal_displ(1:3:length(nodal_displ));
nodal_delta_y = nodal_displ(2:3:length(nodal_displ));
nodal_delta_z = nodal_displ(3:3:length(nodal_displ));
new_nodal_coord = [nodal_delta_x nodal_delta_y nodal_delta_z] + nodal_coordinate;
draw3DMesh(new_nodal_coord, faces);

counter_1 = 3;
displ_mesh = zeros(slab_mesh_meta_data(2)+1, slab_mesh_meta_data(1)+1);

for ii = 1:slab_mesh_meta_data(2)+ 1
    for jj = 1:slab_mesh_meta_data(1)+1
        displ_mesh(ii, jj) = nodal_displ(3*(slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1)*(slab_mesh_meta_data(3)) + counter_1);
        counter_1 = counter_1 + 3;
    end
end
no_nodes_slab = ((slab_mesh_meta_data(1)+1)* (slab_mesh_meta_data(2)+1)*(slab_mesh_meta_data(3)+1));
distinct_x = unique(nodal_coordinate(1:no_nodes_slab, 1), 'rows');
distinct_y = unique(nodal_coordinate(1:no_nodes_slab, 2), 'rows');

figure;
contourf(distinct_x, distinct_y, displ_mesh);
colorbar;

figure;
surf(distinct_x, distinct_y, displ_mesh);
colorbar;
toc

%%Stress
element = (slab_mesh_meta_data(3)-1)*slab_mesh_meta_data(2)*slab_mesh_meta_data(1) ;

distinct_x = unique(nodal_coordinate(1:no_nodes_slab, 1), 'rows');
distinct_y = unique(nodal_coordinate(1:no_nodes_slab, 2), 'rows');
temp_counter = 1;
stress_counter = 1;
stress_mesh_xz = zeros(slab_mesh_meta_data(2)+1, slab_mesh_meta_data(1)+1);

% stress_coord = [slab_dimension(1)/2 0 slab_dimension(3)+col_dimension(1,3) ; 0 slab_dimension(2)/2 slab_dimension(3)+col_dimension(1,3) ; col_top_postion(1,1) col_top_postion(1,2) slab_dimension(3)+col_dimension(1,3) ];
stress_coord = [slab_dimension(1)/2 0  ; 0 slab_dimension(2)/2 ; col_top_postion(1,1) col_top_postion(1,2)];
slab_coord = [ nodal_coordinate(1:(slab_mesh_meta_data(1)+1)*(slab_mesh_meta_data(2)+1),1:2)];
for ii =1:length(stress_coord)
distances(ii,:) = sqrt(sum(bsxfun(@minus, slab_coord, stress_coord(ii,:)).^2,2));
[mi, I] = min(distances(ii,:));
closest(ii,:) = slab_coord(I,:);
% closest(ii,:) = slab_coord(find(distances(ii,:)==min(distances(ii,:))));
end
stress_disp = zeros(3,1);

for ii = 1:slab_mesh_meta_data(2)+ 1
    for jj = 1:slab_mesh_meta_data(1)+1
        if(ii == slab_mesh_meta_data(2) && temp_counter == 0)
           temp_counter = stress_counter;
        end
        if(temp_counter ~= 0 && ii == slab_mesh_meta_data(2)+1)
            stress_counter = temp_counter;
            if(jj == slab_mesh_meta_data(1) + 1)
                stress_mesh_xz(ii, jj) = each_ele_stress(7, 3, element+stress_counter-1);
            else
                stress_mesh_xz(ii, jj) = each_ele_stress(8, 3, element+stress_counter);
                stress_counter = stress_counter + 1;
            end
            temp_counter = stress_counter;
        elseif(jj == slab_mesh_meta_data(2) + 1)
            stress_mesh_xz(ii, jj) = each_ele_stress(5, 3,element+stress_counter-1);
        else
            stress_mesh_xz(ii, jj) = each_ele_stress(6, 3,element+stress_counter);
            stress_counter = stress_counter + 1;
        end
        % Conserdering the three points of stress
        if ((nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),1) == closest(1,1))&(nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),2) == closest(1,2)) )
            stress_disp(1) = stress_mesh_xz(ii, jj);
        elseif ((nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),1) == closest(2,1))&(nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),2) == closest(2,2)) ) 
            stress_disp(2) = stress_mesh_xz(ii, jj);
        elseif ((nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),1) == closest(3,1))&(nodal_coordinate((jj + (slab_mesh_meta_data(1)+ 1)*(ii-1)),2) == closest(3,2)) )
            stress_disp(3) = stress_mesh_xz(ii, jj);
        end
    end
end
figure;
contourf(distinct_x, distinct_y, stress_mesh_xz);
colorbar;

%%strain
temp_counter = 1;
strain_counter = 1;
strain_mesh_xz = zeros(slab_mesh_meta_data(2)+1, slab_mesh_meta_data(1)+1);

for ii = 1:slab_mesh_meta_data(2)+ 1
    for jj = 1:slab_mesh_meta_data(1)+1
        if(ii == slab_mesh_meta_data(2) && temp_counter == 0)
           temp_counter = strain_counter;
        end
        if(temp_counter ~= 0 && ii == slab_mesh_meta_data(2)+1)
            strain_counter = temp_counter;
            if(jj == slab_mesh_meta_data(1) + 1)
                strain_mesh_xz(ii, jj) = each_ele_strain(7,3, element+strain_counter-1);
            else
                strain_mesh_xz(ii, jj) = each_ele_strain(8, 3, element+strain_counter);
                strain_counter = strain_counter + 1;
            end
            temp_counter = strain_counter;
        elseif(jj == slab_mesh_meta_data(2) + 1)
            strain_mesh_xz(ii, jj) = each_ele_strain(5, 3,element+strain_counter-1);
        else
            strain_mesh_xz(ii, jj) = each_ele_strain(6, 3,element+strain_counter);
            strain_counter = strain_counter + 1;
        end
    end
end
figure;
contourf(distinct_x, distinct_y, strain_mesh_xz);
colorbar;