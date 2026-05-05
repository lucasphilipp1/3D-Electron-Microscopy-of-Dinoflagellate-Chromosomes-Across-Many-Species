%make schematic where front half and back half monomers are colored
%differently for LH/RH helices and for a Cholesteric Liquid Crystal (CLC)
%chromosome model

% Parameters for the helix
radius = 2.5;   % Radius of the helix
pitch = -1.5;  % Distance between turns, negative -> left handed helix, positive -> right handed helix
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

total_points = size(xyz_coordinates,1);

%front half and back half of helix are different colors
front_back_half_colormap = [];
for i = 1:1:size(xyz_coordinates(1:1:end,2))
    temp=xyz_coordinates(1:1:end,2);
    if temp(i)>0
        %front_back_half_colormap = [front_back_half_colormap; [0.3293, 0.1227, 0.3707]];	%light purple
        front_back_half_colormap = [front_back_half_colormap; [0.3, 0.3, 0.3]];
    else
        %front_back_half_colormap = [front_back_half_colormap; [0.8884, 0.3309, 1.0000]]; %dark purple
        front_back_half_colormap = [front_back_half_colormap; [0.8, 0.8, 0.8]]; %dark purple
    end
end

figure
axis off
screen = get(0, 'Screensize');
screen(3)=screen(3)/1.75;
set(gcf, 'Position', screen);
hold on
for i = 1:total_points-1
    plot3([xyz_coordinates(i,1), xyz_coordinates(i+1,1)], [xyz_coordinates(i,2), xyz_coordinates(i+1,2)], [xyz_coordinates(i,3), xyz_coordinates(i+1,3)], 'LineWidth', 2, 'Color', [0 0 0]);
end
scatter3(xyz_coordinates(1:1:end,1),xyz_coordinates(1:1:end,2),xyz_coordinates(1:1:end,3),2000,front_back_half_colormap,'filled')
h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
        'MarkerFaceColor', [.8, .8, .8], 'MarkerEdgeColor', 'none');
    %'MarkerFaceColor', [0.8884, 0.3309, 1.0000], 'MarkerEdgeColor', 'none');
h2 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
        'MarkerFaceColor', [.3, .3, .3], 'MarkerEdgeColor', 'none');
    %'MarkerFaceColor', [0.3293, 0.1227, 0.3707], 'MarkerEdgeColor', 'none');
%lgd=legend([h1, h2], {'Front Half', 'Back Half'}, 'Location', 'northeastoutside');
%lgd.FontSize = 32;
%lgd.Box = 'off';
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'Zticklabel',[])
title('Left Handed Helix', FontSize=24);
axis equal
view(0,0) %side view
%view(30,15) %perspective view
%ylim([min(xyz_coordinates(:,2))*3 -0.01]) %show front half only
%ylim([0.01 max(xyz_coordinates(:,2))*3]) %show back half only

% opposite handed helix (right handed)
pitch = 1.5;  % Distance between turns, negative -> left handed helix, positive -> right handed helix

% Create a vector of points along the z-axis
z = linspace(0, num_turns * pitch, total_points);

vertices = [x' y' z'];

% Apply the rotation to the vertices
rotated_vertices = (R * vertices')';

%use for helix plot with tangent vectors
xyz_coordinates = rotated_vertices;

total_points = size(xyz_coordinates,1);

%front half and back half of helix are different colors
front_back_half_colormap = [];
for i = 1:1:size(xyz_coordinates(1:1:end,2))
    temp=xyz_coordinates(1:1:end,2);
    if temp(i)>0
        front_back_half_colormap = [front_back_half_colormap; [0.3293, 0.1227, 0.3707]];	%light purple
    else
        front_back_half_colormap = [front_back_half_colormap; [0.8884, 0.3309, 1.0000]]; %dark purple
    end
end

figure
axis off
screen = get(0, 'Screensize');
screen(3)=screen(3)/1.75;
set(gcf, 'Position', screen);
hold on
for i = 1:total_points-1
    plot3([xyz_coordinates(i,1), xyz_coordinates(i+1,1)], [xyz_coordinates(i,2), xyz_coordinates(i+1,2)], [xyz_coordinates(i,3), xyz_coordinates(i+1,3)], 'LineWidth', 2, 'Color', [0 0 0]);
end
scatter3(xyz_coordinates(1:1:end,1),xyz_coordinates(1:1:end,2),xyz_coordinates(1:1:end,3),2000,front_back_half_colormap,'filled')
h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.8884, 0.3309, 1.0000], 'MarkerEdgeColor', 'none');
h2 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.3293, 0.1227, 0.3707], 'MarkerEdgeColor', 'none');
lgd=legend([h1, h2], {'Front Half', 'Back Half'}, 'Location', 'northeastoutside');
lgd.FontSize = 32;
lgd.Box = 'off';
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'Zticklabel',[])
title('Right Handed Helix', FontSize=24);
axis equal
view(0,0) %side view
%view(30,15) %perspective view
%ylim([min(xyz_coordinates(:,2))*3 -0.01]) %show front half only
%ylim([0.01 max(xyz_coordinates(:,2))*3]) %show back half only

%load CLC model from iter_Cholesteric_HiC.m
chromosome_w_inter_and_intra_disc_loops = readmatrix('chromosome_w_inter_and_intra_disc_loops.csv');
xyz_coordinates=chromosome_w_inter_and_intra_disc_loops;
total_points = size(chromosome_w_inter_and_intra_disc_loops,1);

%front half and back half of CLC are different colors
front_back_half_colormap = [];
for i = 1:1:size(xyz_coordinates(1:1:end,2))
    temp=xyz_coordinates(1:1:end,2);
    if temp(i)>0
        front_back_half_colormap = [front_back_half_colormap; [0.3293, 0.1227, 0.3707]];	%light purple
    else
        front_back_half_colormap = [front_back_half_colormap; [0.8884, 0.3309, 1.0000]]; %dark purple
    end
end

figure
axis off
screen = get(0, 'Screensize');
screen(3)=screen(3)/1.75;
set(gcf, 'Position', screen);
hold on
for i = 1:total_points-1
    plot3([xyz_coordinates(i,1), xyz_coordinates(i+1,1)], [xyz_coordinates(i,2), xyz_coordinates(i+1,2)], [xyz_coordinates(i,3), xyz_coordinates(i+1,3)], 'LineWidth', 2, 'Color', [0 0 0]);
end
scatter3(xyz_coordinates(1:1:end,1),xyz_coordinates(1:1:end,2),xyz_coordinates(1:1:end,3),400,front_back_half_colormap,'filled')
% Plot dummy points for color legend entries
h1 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.8884, 0.3309, 1.0000], 'MarkerEdgeColor', 'none');
h2 = plot3(nan, nan, nan, 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', [0.3293, 0.1227, 0.3707], 'MarkerEdgeColor', 'none');
lgd=legend([h1, h2], {'Front Half', 'Back Half'}, 'Location', 'northeast');
lgd.FontSize = 32;
lgd.Box = 'off';
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
set(gca,'Zticklabel',[])
axis equal
%view(0,0) %side view
view(30,15) %perspective view
%ylim([min(xyz_coordinates(:,2))*3 -0.01]) %show front half only
%ylim([0.01 max(xyz_coordinates(:,2))*3]) %show back half only

