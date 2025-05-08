clc
clear

% Parameters for the helix
radius = 10;   % Radius of the helix
pitch = -10;  % Distance between turns, negative -> left handed helix, positive -> right handed helix
num_turns = 5; % Number of turns
num_points_per_turn = 40; % Number of points per turn

% Calculate the total number of points
total_points = num_turns * num_points_per_turn;

% Create a vector of points along the z-axis
z = linspace(0, num_turns * pitch, total_points);

% Create a vector of angles for the helix
theta = linspace(0, num_turns * 2 * pi, total_points);

% Calculate the x and y coordinates of the helix
x = radius * cos(theta);
y = radius * sin(theta);

vertices = [x' y' z'];

% Generate random angles for rotation around x, y, and z axes
theta_x = rand * 2 * pi;
theta_y = rand * 2 * pi;
theta_z = rand * 2 * pi;

% or no rotation
theta_x = 0;
theta_y = 0;
theta_z = 0;

% Create rotation matrices for each axis
Rx = [
    1 0 0;
    0 cos(theta_x) -sin(theta_x);
    0 sin(theta_x) cos(theta_x)
    ];

Ry = [
    cos(theta_y) 0 sin(theta_y);
    0 1 0;
    -sin(theta_y) 0 cos(theta_y)
    ];

Rz = [
    cos(theta_z) -sin(theta_z) 0;
    sin(theta_z) cos(theta_z) 0;
    0 0 1
    ];

% Combine the rotation matrices
R = Rx * Ry * Rz;

% Apply the rotation to the vertices
rotated_vertices = (R * vertices')';

%use for helix plot with tangent vectors
xyz_coordinates = rotated_vertices;

%use for binary image stack generation
%xyz_coordinates=round(rotated_vertices)-(round(min(min(rotated_vertices)))-1);

monomer_index=1:1:size(xyz_coordinates,1);
query_index=1:0.1:size(xyz_coordinates,1);

xx = interp1(monomer_index,xyz_coordinates(:,1),query_index,'spline');
yy = interp1(monomer_index,xyz_coordinates(:,2),query_index,'spline');
zz = interp1(monomer_index,xyz_coordinates(:,3),query_index,'spline');

bond_vectors = [];

for j = 2:1:size(xyz_coordinates,1)-1
    left_point_x=xx(find(query_index==monomer_index(j))-1);
    right_point_x=xx(find(query_index==monomer_index(j))+1);

    left_point_y=yy(find(query_index==monomer_index(j))-1);
    right_point_y=yy(find(query_index==monomer_index(j))+1);

    left_point_z=zz(find(query_index==monomer_index(j))-1);
    right_point_z=zz(find(query_index==monomer_index(j))+1);

    tx=right_point_x-left_point_x;
    ty=right_point_y-left_point_y;
    tz=right_point_z-left_point_z;

    bond_vectors = [bond_vectors; [tx./norm([tx ty tz]) ty./norm([tx ty tz]) tz./norm([tx ty tz])]];
end

bond_vectors = [[0 0 0]; bond_vectors; [0 0 0]];

resolution = 5000;
chromosome = [xx' yy' zz'];
chromosome = flipud(chromosome);

numPoints = size(chromosome,1);

MyColor = linspace(1,numPoints,numPoints)';
% create a connectivity matrix
Faces = [1:(numPoints-1); 2:numPoints]';

% Plot the original helix
figure
screen = get(0, 'Screensize');
screen(3)=screen(3)/1.75;
set(gcf, 'Position', screen);
% Apply interpolation for each x,y and z
hold on
%plot3(xx,yy,zz,'Color', [0 0 0], LineWidth=15)
c = colormap('jet');

idx=1:10:size(xx,2);
color_idx=round(linspace(1,256,100));
color_idx = fliplr(color_idx);

for i = 1:total_points-1
    % Interpolate between two points (x(i), y(i), z(i)) and (x(i+1), y(i+1), z(i+1))
    % Plot a line segment between successive points with the appropriate color
    plot3([xx(idx(i)), xx(idx(i+1))], [yy(idx(i)), yy(idx(i+1))], [zz(idx(i)), zz(idx(i+1))], 'LineWidth', 2, 'Color', c(color_idx(i),:)); 
end

scatter3(xx(1:10:end),yy(1:10:end),zz(1:10:end),1000,z,'filled')
colormap jet
%start=300;
%stop=370;
%quiver3(xyz_coordinates(start/10,1)-bond_vectors(start/10,1),xyz_coordinates(start/10,2)-bond_vectors(start/10,2),xyz_coordinates(start/10,3)-bond_vectors(start/10,3),bond_vectors(start/10,1),bond_vectors(start/10,2),bond_vectors(start/10,3),12,'r',"LineWidth", 5, "MaxHeadSize",45)
%quiver3(xyz_coordinates(stop/10,1)-bond_vectors(stop/10,1),xyz_coordinates(stop/10,2)-bond_vectors(stop/10,2),xyz_coordinates(stop/10,3)-bond_vectors(stop/10,3),bond_vectors(stop/10,1),bond_vectors(stop/10,2),bond_vectors(stop/10,3),12,'r',"LineWidth", 5, "MaxHeadSize",45)
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'Zticklabel',[])
axis equal
view(90,0)

colormap jet
axis equal
% caxis([min(MyColor) max(MyColor)])
% c = colorbar;
% c.Position = c.Position - [.1 0 0 0];
% c.Ticks = linspace(0, size(chromosome,1), round(size(chromosome,1)/500)+1);
% c.TickLabels = num2cell(linspace(0, 2*10^7, round(size(chromosome,1)/500)+1));
% c.Label.String = 'primary sequence [bp]';
% c.FontSize = 32;
% patch('Faces', Faces(:,:) ,'Vertices', chromosome(:,:) ,'FaceColor', 'none', 'FaceVertexCData', MyColor(:,:) ,'EdgeColor','interp' ,'LineWidth',25, 'FaceAlpha',.5,'EdgeAlpha',.5);

%view(45,22.5)
view(0,0)
% 
% xlim([min(chromosome(:,1))*1.6 max(chromosome(:,1))*1.6])
% ylim([min(chromosome(:,2))*1.6 max(chromosome(:,2))*1.6])
% %ylim([0 max(yy)*1.6])
% ylim([min(yy)*1.6 0])
% zlim([min(chromosome(:,3))*1.1 max(chromosome(:,3))*1.1])
% 
% % Determine the size of the image stack based on the maximum X, Y, and Z coordinates
% max_x = max(max(xyz_coordinates));
% max_y = max(max(xyz_coordinates));
% max_z = max(max(xyz_coordinates));

% % Initialize the binary image stack
% binary_stack = false(max_y, max_x, max_z);

% % Convert XYZ coordinates to binary image stack
% for i = 1:size(xyz_coordinates, 1)
%     x = xyz_coordinates(i, 1);
%     y = xyz_coordinates(i, 2);
%     z = xyz_coordinates(i, 3);
%     binary_stack(y, x, z) = true;
% end
%
% % Display one of the slices (optional)
% for i = 1:1:size(binary_stack,3)
% imwrite(binary_stack(:, :, i),sprintf("/Users/lucasphilipp/Downloads/helix/helix_%i.png",i));
% end
