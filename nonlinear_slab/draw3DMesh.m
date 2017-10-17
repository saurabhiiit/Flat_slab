function draw3DMesh(nodal_coordinate, faces)
%**************************************************************************
% Uses patch to draw 3D mesh for the given wall.
%**************************************************************************
%
% Input parameters:
% nodal_coordinate - Nodal coordinate matrix.
% faces            - Face matrix containing face connection for each
%                    elements together for plotting using patch.

figure;
patch('Vertices',nodal_coordinate,'Faces',faces,...
    'FaceColor',[0.7 0.7 0.7]);
 axis off;
% title('3D model of flat slab','FontSize', 30)
% legend('EW')
set(gcf,'PaperPosition',[0 0 10 7])
% print -djpeg -f4001 -r300
 % axis equal;
% cameratoolbar('SetCoordSys','z');
% cameratoolbar('setmode','orbit');
% ax = gca;
% ax.CameraPosition = [-179534.00986074534 128291.59093640426 128119.23378678535];
% ax.CameraUpVector = [0.40957602214449595 -0.28678821817552286 0.8660254037844387];
% ax.CameraViewAngle = 1.4342442304781926;
end