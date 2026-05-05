clc
clear

folderPath = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/all crescents final';
files = dir(folderPath);
files = files(~[files.isdir]);
fileNames = fullfile(folderPath, {files.name})';
fileNames = fileNames(~contains(fileNames, '.DS_Store'));
numFiles = size(fileNames,1);

resolution = 4; %voxel size in nanometers

crescent_arc_lengths = [];
crescent_widths = [];

for l = 1:1:numFiles
    % Load 3D binary TIFF
    info = imfinfo(fileNames{l});
    num_slices = numel(info);
    volume = false(info(1).Height, info(1).Width, num_slices);
    for k = 1:num_slices
        volume(:,:,k) = imread(fileNames{l}, 'Index', k) > 0;
    end

    volume = imfill(volume,'holes');

    % 3D skeletonize
    skel3D = bwskel(volume, 'MinBranchLength', 10); %MinBranchLength needs to be tuned

    [y,x,z] = ind2sub(size(skel3D), find(skel3D));
    pts = [x, y, z];
    numPts = size(pts,1);

    G = graph;
    [idx, dist] = rangesearch(pts, pts, 3); %the search range needs to be tuned
    for i = 1:numPts
        neighbors = idx{i};
        separation = dist{i};
        for j = 1:length(neighbors)
            n = neighbors(j);
            if n > i
                G = addedge(G, i, n, separation(j));
            end
        end
    end

    n = numnodes(G);
    D = distances(G);  % compute shortest path between all node pairs
    D(isinf(D)) = -inf;

    % Find max shortest path length and indices
    [max_dist, linear_idx] = max(D(:));
    [start_node, end_node] = ind2sub([n n], linear_idx);

    % Extract the longest shortest path
    path = shortestpath(G, start_node, end_node);

    % Compute arclength
    total_length = 0;
    for i = 1:length(path)-1
        total_length = total_length + norm(pts(path(i+1), :) - pts(path(i), :));
    end
    crescent_arc_lengths = [crescent_arc_lengths; total_length*resolution];
    %fprintf('3D arc length: %.2f voxels\n', total_length);

    % Compute width
    D3 = bwdist(~volume);
    widths = zeros(length(path),1);
    for i = 1:length(path)
        coords = pts(path(i), :);
        widths(i) = 2 * D3(coords(2), coords(1), coords(3));
    end
    crescent_mean_width = mean(widths);
    crescent_widths = [crescent_widths; crescent_mean_width*resolution];
    %fprintf('Mean width: %.2f voxels\n', mean_width);

    %Visualize width along centerline
    % figure;
    % hold on
    % scatter3(pts(path,1), pts(path,2), pts(path,3), 40, widths, 'filled');
    % scatter3(pts(start_node,1), pts(start_node,2), pts(start_node,3), 100, 'r', 'filled');
    % scatter3(pts(end_node,1), pts(end_node,2), pts(end_node,3), 100, 'g', 'filled');
    % 
    % fv = isosurface(volume, 0.5);
    % p = patch(fv);
    % p.FaceColor = [0.8 0.8 0.8];  % Light gray surface
    % p.EdgeColor = 'none';
    % p.FaceAlpha = 0.2;
    % 
    % caxis([0 30])
    % 
    % cb=colorbar;
    % cb.Label.String = ('Width [nm]');
    % cb.Label.FontSize = 20;
    % cb.FontSize = 20;
    % xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]');
    % 
    % ax = gca;
    % ax.FontSize = 20;
    % 
    % colormap(jet);
    % 
    % grid on;
    % view([-0.2 1 0.1]);
    % axis vis3d;
    % hold off
    % [~, filename, ~] = fileparts(fileNames{l});
    % print(gcf, [filename, '.png'], '-dpng', '-r300');

    l
end

folderPath = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/all rods final';
files = dir(folderPath);
files = files(~[files.isdir]);
fileNames = fullfile(folderPath, {files.name})';
fileNames = fileNames(~contains(fileNames, '.DS_Store'));
numFiles = size(fileNames,1);

rod_lengths = [];
rod_widths = [];

for l = 1:1:numFiles
    % Load 3D binary TIFF
    info = imfinfo(fileNames{l});
    num_slices = numel(info);
    volume = false(info(1).Height, info(1).Width, num_slices);
    for k = 1:num_slices
        volume(:,:,k) = imread(fileNames{l}, 'Index', k) > 0;
    end

    volume = imfill(volume,'holes');

    % 3D skeletonize
    skel3D = bwskel(volume, 'MinBranchLength', 10); %MinBranchLength needs to be tuned

    [y,x,z] = ind2sub(size(skel3D), find(skel3D));
    pts = [x, y, z];
    numPts = size(pts,1);

    G = graph;
    [idx, dist] = rangesearch(pts, pts, 3); %the search range needs to be tuned
    for i = 1:numPts
        neighbors = idx{i};
        separation = dist{i};
        for j = 1:length(neighbors)
            n = neighbors(j);
            if n > i
                G = addedge(G, i, n, separation(j));
            end
        end
    end

    n = numnodes(G);
    D = distances(G);  % compute shortest path between all node pairs
    D(isinf(D)) = -inf;

    % Find max shortest path length and indices
    [max_dist, linear_idx] = max(D(:));
    [start_node, end_node] = ind2sub([n n], linear_idx);

    % Extract the longest shortest path
    path = shortestpath(G, start_node, end_node);

    % Compute length
    total_length = 0;
    for i = 1:length(path)-1
        total_length = total_length + norm(pts(path(i+1), :) - pts(path(i), :));
    end
    rod_lengths = [rod_lengths; total_length*resolution];
    %fprintf('3D arc length: %.2f voxels\n', total_length);

    % Compute width
    D3 = bwdist(~volume);
    widths = zeros(length(path),1);
    for i = 1:length(path)
        coords = pts(path(i), :);
        widths(i) = 2 * D3(coords(2), coords(1), coords(3));
    end
    rod_mean_width = mean(widths);
    rod_widths = [rod_widths; rod_mean_width*resolution];
    %fprintf('Mean width: %.2f voxels\n', mean_width);

    l
end

folderPath = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/all toroids final'; %folder path with all the 3D binary tiff crescent volumes
files = dir(folderPath);
files = files(~[files.isdir]);
fileNames = fullfile(folderPath, {files.name})';
fileNames = fileNames(~contains(fileNames, '.DS_Store'));
numFiles = size(fileNames,1);

toroid_circumference = [];
toroid_widths = [];

for l = 1:1:numFiles
    % Load 3D binary TIFF
    info = imfinfo(fileNames{l});
    num_slices = numel(info);
    volume = false(info(1).Height, info(1).Width, num_slices);
    for k = 1:num_slices
        volume(:,:,k) = imread(fileNames{l}, 'Index', k) > 0;
    end

    volume = imfill(volume,'holes');

    %toroid_vol = [toroid_vol; sum(I,"all")*4^3]; %in nm^3
    s = regionprops3(volume,"PrincipalAxisLength");
    toroid_widths = [toroid_widths; 4*s.PrincipalAxisLength(3)]; %in nm
    toroid_circumference = [toroid_circumference; (4*(s.PrincipalAxisLength(1)+s.PrincipalAxisLength(1))/2)*pi]; %in nm

    l
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chromosome_vol_Brandtodinium_nutricula = [];
chromosome_vol_Symbiodinium_pilosum = [];
chromosome_vol_Symbiodinium_microadriaticum = [];
chromosome_vol_Symbiodinium_minutum = [];
chromosome_vol_Crypthecodinium_cohnii = [];
chromosome_vol_Ensiculifera_tyrrhenica = [];
chromosome_vol_Ross_Sea_Dinoflagellate = [];


chromosome_dimensions = [];

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Brandtodinium nutricula/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        if ii == 1
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*4^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 2
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*4^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 3
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*4^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        end
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium microadriaticum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_microadriaticum = [chromosome_vol_Symbiodinium_microadriaticum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium minutum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_minutum = [chromosome_vol_Symbiodinium_minutum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium pilosum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_pilosum = [chromosome_vol_Symbiodinium_pilosum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Crypthecodinium cohnii/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Crypthecodinium_cohnii = [chromosome_vol_Crypthecodinium_cohnii; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Ensiculifera_tyrrhenica/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Ensiculifera_tyrrhenica = [chromosome_vol_Ensiculifera_tyrrhenica; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Ross Sea Dinoflagellate/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Ross_Sea_Dinoflagellate = [chromosome_vol_Ross_Sea_Dinoflagellate; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        chromosome_dimensions = [chromosome_dimensions; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
edges = linspace(0, 3000, 17);
edges2 = linspace(0, 250, 11);

figure('Units', 'normalized', 'Position', [0.2, 0.1, 0.3, 0.7])
subplot(2,1,1)
pos1 = get(gca, 'Position');
hold on

%% Crescent Width and Arc Length
[N1, e1] = histcounts(crescent_widths, edges2);
[N2, e2] = histcounts(crescent_arc_lengths, edges);
e1 = e1(2:end) - diff(e1(1:2))/2;
e2 = e2(2:end) - diff(e2(1:2))/2;
N1(N1 == 0) = -1;
N2(N2 == 0) = -1;
%plot(e1, N1, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1)
plot(e2, N2, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
%text(e1(1), N1(1)-5, 'c', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
for i = 1:length(e2)
    if N2(i) > 0
        text(e2(i), N2(i), 'c', 'FontSize', 24, 'Color', [0.75 0.75 0.75], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

%% Toroid Width and Circumference
[N3, e3] = histcounts(toroid_widths, edges2);
[N4, e4] = histcounts(toroid_circumference, edges);
e3 = e3(2:end) - diff(e3(1:2))/2;
e4 = e4(2:end) - diff(e4(1:2))/2;
N3(N3 == 0) = -1;
N4(N4 == 0) = -1;
%plot(e3, N3, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1)
plot(e4, N4, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1)
%text(e3(1), N3(1)+2, 'o', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
for i = 1:length(e4)
    if N4(i) > 0
        text(e4(i), N4(i), 'o', 'FontSize', 24, 'Color', [0.75 0.75 0.75], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

%% Rod Width and Length
[N5, e5] = histcounts(rod_widths, edges2);
[N6, e6] = histcounts(rod_lengths, edges);
e5 = e5(2:end) - diff(e5(1:2))/2;
e6 = e6(2:end) - diff(e6(1:2))/2;
N5(N5 == 0) = -1;
N6(N6 == 0) = -1;
%plot(e5, N5, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.5)
plot(e6, N6, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5)
%text(e5(1), N5(1), 'r', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
for i = 1:length(e6)
    if N6(i) > 0
        text(e6(i), N6(i), 'r', 'FontSize', 24, 'Color', [0.75 0.75 0.75], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

% % Legend handles
% h_rod_dark = plot(NaN, NaN,'-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1);
% h_crescent_dark = plot(NaN, NaN,'-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1);
% h_toroid_dark = plot(NaN, NaN,'-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1);
% 
% h_rod_light = plot(NaN, NaN,'-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);    
% h_crescent_light = plot(NaN, NaN,'-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
% h_toroid_light = plot(NaN, NaN,'-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
% 
% legend([h_rod_dark, h_crescent_dark, h_toroid_dark, ...
%         h_rod_light, h_crescent_light, h_toroid_light], ...
%        {'Rod Width', 'Crescent Width', 'Toroid Width', ...
%         'Rod Length', 'Crescent Arc Length', 'Toroid Circumference'}, ...
%        'Location', 'northeast');
% legend boxoff

ylabel('Count', 'FontSize', 22.5)
xlim([0 3000])
ylim([0 120])
ylim([0 100])

ax = gca;
ax.FontSize = 22.5;

%text(1375, 111, 'r', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center');
%text(1375, 100, 'c', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center');
%text(1375, 89, 'o', 'FontSize', 24, 'Color', [0.25 0.25 0.25], 'HorizontalAlignment', 'center');
%text(1375, 78, 'r', 'FontSize', 24, 'Color', [0.75 0.75 0.75], 'HorizontalAlignment', 'center');
%text(1375, 67, 'c', 'FontSize', 24, 'Color', [0.75 0.75 0.75], 'HorizontalAlignment', 'center');
%text(1375, 56, 'o', 'FontSize', 24, 'Color', [0.75 0.75 0.75], 'HorizontalAlignment', 'center');

%% === Subplot 2: Chromosome Dimensions ===
subplot(2,1,2)
pos2 = get(gca, 'Position');

hold on

[N7, e7] = histcounts((chromosome_dimensions(:,2)+chromosome_dimensions(:,3))/2, edges);
[N8, e8] = histcounts(chromosome_dimensions(:,1), edges);
e7 = e7(2:end) - diff(e7(1:2))/2;
e8 = e8(2:end) - diff(e8(1:2))/2;
N7(N7 == 0) = -1;
N8(N8 == 0) = -1;

plot(e7, N7, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 3)
plot(e8, N8, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 3)

% Axis settings
xlabel('DNA Object Dimensions [nm]', 'FontSize', 22.5)
ylabel('Count', 'FontSize', 22.5)
xlim([0 3000])
ylim([0 800])
ax = gca;
ax.FontSize = 22.5;

% Chromosome legend
h_chr_dark = plot(NaN, NaN,'-', 'Color', [0.25 0.25 0.25], 'LineWidth', 3);
h_chr_light = plot(NaN, NaN,'-', 'Color', [0.75 0.75 0.75], 'LineWidth', 3);
legend([h_chr_dark, h_chr_light], {'Chromosome Width','Chromosome Length'}, 'Location', 'northeast');
legend boxoff

% Bring plots closer together vertically
pos1(2) = pos1(2) - 0.075/2;
pos2(4) = pos2(4) + 0.075/2;

pos1(1) = pos1(1) + 0.075/2;
pos2(1) = pos2(1) + 0.075/2;

set(gca, 'Position', pos2)  % Update subplot 2 position

subplot(2,1,1)
set(gca, 'Position', pos1)  % Update subplot 1 position

%%%

figure
hold on
plot(e1, N1, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 3)
plot(e3, N3, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 3)
plot(e5, N5, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 3)
for i = 1:length(e1)
    if N1(i) > 0
        text(e1(i), N1(i), 'c', 'FontSize', 40, 'Color', [0.25 0.25 0.25], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
for i = 1:length(e3)
    if N3(i) > 0
        text(e3(i), N3(i), 'o', 'FontSize', 40, 'Color', [0.25 0.25 0.25], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
for i = 1:length(e5)
    if N5(i) > 0
        text(e5(i), N5(i), 'r', 'FontSize', 40, 'Color', [0.25 0.25 0.25], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
ax = gca;
ax.FontSize = 32;
xlim([0 250])
ylim([0 40])