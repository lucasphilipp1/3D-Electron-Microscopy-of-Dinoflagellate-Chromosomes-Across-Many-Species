clc
clear

% Parameters for the helix
radius = 5;   % Radius of the helix
pitch = 10;  % Distance between turns
num_turns = 5; % Number of turns
num_points_per_turn = 100; % Number of points per turn

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

%%Your points
xyz_coordinates=round(rotated_vertices)-(round(min(min(rotated_vertices)))-1);

% Plot the original helix
figure
monomer_index=1:1:size(xyz_coordinates,1);
query_index=1:10:size(xyz_coordinates,1);
% Apply interpolation for each x,y and z
xx = interp1(monomer_index,xyz_coordinates(:,1),query_index,'spline');
yy = interp1(monomer_index,xyz_coordinates(:,2),query_index,'spline');
zz = interp1(monomer_index,xyz_coordinates(:,3),query_index,'spline');
plot3(xx,yy,zz,'-k', LineWidth=2)
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[]) 
set(gca,'Zticklabel',[])

% Determine the size of the image stack based on the maximum X, Y, and Z coordinates
max_x = max(max(xyz_coordinates));
max_y = max(max(xyz_coordinates));
max_z = max(max(xyz_coordinates));

% Initialize the binary image stack
binary_stack = false(max_y, max_x, max_z);

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
