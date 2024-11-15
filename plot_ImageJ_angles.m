clc
clear

temp = [];
sto = [];

s_pilosum_angles = [];
s_pilosum_vol_no_angles = [];

b_nutricula_angles = [];
b_nutricula_vol_no_angles = [];

s_minutum_angles = [];
s_minutum_vol_no_angles = [];

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP cell 1.csv");
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_pilosum_cell_1_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_pilosum_angles(:,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/',file_names{i}));
I=I./max(I);
s_pilosum_vol_no_angles = [s_pilosum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP cell 2.csv");
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_pilosum_cell_2_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_pilosum_angles(s_pilosum_cell_1_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/',file_names{i}));
I=I./max(I);
s_pilosum_vol_no_angles = [s_pilosum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/SP cell 3.csv");
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_pilosum_cell_3_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_pilosum_angles(s_pilosum_cell_1_chr_annotated+s_pilosum_cell_2_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/',file_names{i}));
I=I./max(I);
s_pilosum_vol_no_angles = [s_pilosum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN cell 1.csv");
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
b_nutricula_cell_1_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*10^3]; %in nm^3
end
b_nutricula_angles(:,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/',file_names{i}));
I=I./max(I);
b_nutricula_vol_no_angles = [b_nutricula_vol_no_angles; sum(I,"all")*10^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN cell 2.csv");
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
b_nutricula_cell_2_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium cell 2 stack chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
b_nutricula_angles(b_nutricula_cell_1_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium cell 2 stack chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium cell 2 stack chromosomes/',file_names{i}));
I=I./max(I);
b_nutricula_vol_no_angles = [b_nutricula_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/BN cell 3.csv");
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
b_nutricula_cell_3_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 3 chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*2^3]; %in nm^3
end
b_nutricula_angles(b_nutricula_cell_1_chr_annotated+b_nutricula_cell_2_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 3 chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 3 chromosomes/',file_names{i}));
I=I./max(I);
b_nutricula_vol_no_angles = [b_nutricula_vol_no_angles; sum(I,"all")*2^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin cell 1.csv");
s_minutum_angles = [s_minutum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_minutum_cell_1_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 1 chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_minutum_angles(:,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 1 chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 1 chromosomes/',file_names{i}));
I=I./max(I);
s_minutum_vol_no_angles = [s_minutum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin cell 2.csv");
s_minutum_angles = [s_minutum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_minutum_cell_2_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 2 chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_minutum_angles(s_minutum_cell_1_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 2 chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 2 chromosomes/',file_names{i}));
I=I./max(I);
s_minutum_vol_no_angles = [s_minutum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%%%

ImageJ_angles = importdata("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/SMin cell 3.csv");
s_minutum_angles = [s_minutum_angles; ImageJ_angles.data(1:2:end,5) ImageJ_angles.data(2:2:end,5) zeros(size(ImageJ_angles.data(1:2:end,5)))]; %back front volume
s_minutum_cell_3_chr_annotated = size(ImageJ_angles.data(1:2:end,5),1);

for i = 2:2:size(ImageJ_angles.textdata(:,2),1)
ROI_with_angle = regexp(ImageJ_angles.textdata(i,2), '\d+', 'match');
sto = [sto; strcat('ROI',{' '},ROI_with_angle{1},'.tiff')];
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 3 chromosomes/','ROI',{' '},ROI_with_angle{1},'.tiff'));
I=I./max(I);
temp = [temp; sum(I,"all")*4^3]; %in nm^3
end
s_minutum_angles(s_minutum_cell_1_chr_annotated+s_minutum_cell_2_chr_annotated+1:end,3) = temp;
temp=[];

files = dir('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 3 chromosomes/*.tiff');
file_names={files.name};
file_names = file_names(~ismember(file_names, sto));
for i = 1:1:size(file_names,2)
I = tiffreadVolume(strcat('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 3 chromosomes/',file_names{i}));
I=I./max(I);
s_minutum_vol_no_angles = [s_minutum_vol_no_angles; sum(I,"all")*4^3]; %in nm^3
end
sto=[];

%mean angle

figure
hold on
histogram((s_pilosum_angles(:,1)+s_pilosum_angles(:,2))./2,15, 'EdgeColor', [0.4940 0.1840 0.5560], 'FaceColor', [0.4940 0.1840 0.5560])
histogram((s_minutum_angles(:,1)+s_minutum_angles(:,2))./2,15, 'EdgeColor', [0.9290 0.6940 0.1250], 'FaceColor', [0.9290 0.6940 0.1250])
histogram((b_nutricula_angles(:,1)+b_nutricula_angles(:,2))./2,15, 'EdgeColor', [0.6350 0.0780 0.1840], 'FaceColor', [0.6350 0.0780 0.1840])
xline(90, '--', '90$^{\circ}$','FontSize', 24,'Interpreter','latex')
xlim([0 180])
xlabel('$<\theta>$ [Degrees]','FontSize', 24, 'Interpreter','latex')
ylabel('Count','FontSize', 24)
lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northwest', 'FontAngle', 'italic');
legend boxoff
ax = gca;
ax.FontSize = 16;

no_angle_extracted=zeros(3,1);
angle_extracted=zeros(3,1);
source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium pilosum/cell 1';
no_angle_extracted(1) = no_angle_extracted(1) + length(dir([source_dir, '/*.tiff']))-s_pilosum_cell_1_chr_annotated;
angle_extracted(1) = angle_extracted(1) + s_pilosum_cell_1_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium pilosum/cell 2';
no_angle_extracted(1) = no_angle_extracted(1) + length(dir([source_dir, '/*.tiff']))-s_pilosum_cell_2_chr_annotated;
angle_extracted(1) = angle_extracted(1) + s_pilosum_cell_2_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium pilosum/cell 3';
no_angle_extracted(1) = no_angle_extracted(1) + length(dir([source_dir, '/*.tiff']))-s_pilosum_cell_3_chr_annotated;
angle_extracted(1) = angle_extracted(1) + s_pilosum_cell_3_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium minutum/cell 1';
no_angle_extracted(2) = no_angle_extracted(2) + length(dir([source_dir, '/*.tiff']))-s_minutum_cell_1_chr_annotated;
angle_extracted(2) = angle_extracted(2) + s_minutum_cell_1_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium minutum/cell 2';
no_angle_extracted(2) = no_angle_extracted(2) + length(dir([source_dir, '/*.tiff']))-s_minutum_cell_2_chr_annotated;
angle_extracted(2) = angle_extracted(2) + s_minutum_cell_2_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium minutum/cell 3';
no_angle_extracted(2) = no_angle_extracted(2) + length(dir([source_dir, '/*.tiff']))-s_minutum_cell_3_chr_annotated;
angle_extracted(2) = angle_extracted(2) + s_minutum_cell_3_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Brandtodinium nutricula/cell 1';
no_angle_extracted(3) = no_angle_extracted(3) + length(dir([source_dir, '/*.tiff']))-b_nutricula_cell_1_chr_annotated;
angle_extracted(3) = angle_extracted(3) + b_nutricula_cell_1_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Brandtodinium nutricula/cell 2';
no_angle_extracted(3) = no_angle_extracted(3) + length(dir([source_dir, '/*.tiff']))-b_nutricula_cell_2_chr_annotated;
angle_extracted(3) = angle_extracted(3) + b_nutricula_cell_2_chr_annotated;

source_dir = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Brandtodinium nutricula/cell 3';
no_angle_extracted(3) = no_angle_extracted(3) + length(dir([source_dir, '/*.tiff']))-b_nutricula_cell_3_chr_annotated;
angle_extracted(3) = angle_extracted(3) + b_nutricula_cell_3_chr_annotated;

figure
b = bar([angle_extracted'; no_angle_extracted']);
b(1).FaceColor = 'flat';
b(1).CData(1,:) = [0.4940 0.1840 0.5560];
b(1).CData(2,:) = [0.4940 0.1840 0.5560];

b(2).FaceColor = 'flat';
b(2).CData(1,:) = [0.9290 0.6940 0.1250];
b(2).CData(2,:) = [0.9290 0.6940 0.1250];

b(3).FaceColor = 'flat';
b(3).CData(1,:) = [0.6350 0.0780 0.1840];
b(3).CData(2,:) = [0.6350 0.0780 0.1840];

set(gca,'XTickLabel',{'Angle Extracted', 'No Angle Extracted'});
%lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northwest', 'FontAngle', 'italic');
%legend boxoff
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;

%angle difference

figure
hold on
histogram(s_pilosum_angles(:,1)-s_pilosum_angles(:,2),15, 'EdgeColor', [0.4940 0.1840 0.5560], 'FaceColor', [0.4940 0.1840 0.5560])
histogram(s_minutum_angles(:,1)-s_minutum_angles(:,2),15, 'EdgeColor', [0.9290 0.6940 0.1250], 'FaceColor', [0.9290 0.6940 0.1250])
histogram(b_nutricula_angles(:,1)-b_nutricula_angles(:,2),15, 'EdgeColor', [0.6350 0.0780 0.1840], 'FaceColor', [0.6350 0.0780 0.1840])
xline(0, '--', '0$^{\circ}$','FontSize', 24,'Interpreter','latex')
xlim([-90 90])
xlabel('$\Delta\theta$ [Degrees]','FontSize', 24, 'Interpreter','latex')
ylabel('Count','FontSize', 24)
lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northwest', 'FontAngle', 'italic');
legend boxoff
ax = gca;
ax.FontSize = 16;

%%%

a = calc_percentile([s_pilosum_angles(:,3); s_pilosum_vol_no_angles],s_pilosum_angles(:,3));
b = calc_percentile([s_minutum_angles(:,3); s_minutum_vol_no_angles],s_minutum_angles(:,3));
c = calc_percentile([b_nutricula_angles(:,3); b_nutricula_vol_no_angles],b_nutricula_angles(:,3));
d = calc_percentile([s_pilosum_angles(:,3); s_pilosum_vol_no_angles],s_pilosum_vol_no_angles);
e = calc_percentile([s_minutum_angles(:,3); s_minutum_vol_no_angles],s_minutum_vol_no_angles);
f = calc_percentile([b_nutricula_angles(:,3); b_nutricula_vol_no_angles],b_nutricula_vol_no_angles);

x = 0:10:100;

a=sort(a);
b=sort(b);
c=sort(c);
d=sort(d);
e=sort(e);
f=sort(f);
d_a=discretize(a,x);
d_b=discretize(b,x);
d_c=discretize(c,x);
d_d=discretize(d,x);
d_e=discretize(e,x);
d_f=discretize(f,x);
N_a = zeros(10,1);
N_b = zeros(10,1);
N_c = zeros(10,1);
N_d = zeros(10,1);
N_e = zeros(10,1);
N_f = zeros(10,1);

for i = 1:1:10
N_a(i) = length(find(d_a==i))/(length(find(d_a==i))+length(find(d_d==i)));
N_d(i) = length(find(d_d==i))/(length(find(d_a==i))+length(find(d_d==i)));

N_b(i) = length(find(d_b==i))/(length(find(d_b==i))+length(find(d_e==i)));
N_e(i) = length(find(d_e==i))/(length(find(d_b==i))+length(find(d_e==i)));

N_c(i) = length(find(d_c==i))/(length(find(d_c==i))+length(find(d_f==i)));
N_f(i) = length(find(d_f==i))/(length(find(d_c==i))+length(find(d_f==i)));
end

figure
hold on
x = x(2:end) - (x(2)-x(1))/2;
plot(x,N_a, 'Color', [0.4940 0.1840 0.5560], LineWidth=2)
plot(x,N_b, 'Color', [0.9290 0.6940 0.1250], LineWidth=2)
plot(x,N_c, 'Color', [0.6350 0.0780 0.1840], LineWidth=2)
plot(x,N_d,'--', 'Color', [0.4940 0.1840 0.5560], LineWidth=2)
plot(x,N_e,'--', 'Color', [0.9290 0.6940 0.1250], LineWidth=2)
plot(x,N_f,'--', 'Color', [0.6350 0.0780 0.1840], LineWidth=2)
xlim([0 100])
ylim([0 1])
xlabel('Chromosome Volume Percentile','FontSize', 24)
ylabel('Probability','FontSize', 24)
lgd=legend({'\textit{S. pilosum} Angle Extracted','\textit{S. minutum} Angle Extracted','\textit{B. nutricula} Angle Extracted','\textit{S. pilosum} No Angle Extracted','\textit{S. minutum} No Angle Extracted','\textit{B. nutricula} No Angle Extracted'},'Interpreter','Latex','Location', 'eastoutside');
legend boxoff
ax = gca;
ax.FontSize = 16;

percentile_vol_pilosum = calc_percentile(s_pilosum_angles(:,3),s_pilosum_angles(:,3));
percentile_vol_minutum = calc_percentile(s_minutum_angles(:,3),s_minutum_angles(:,3));
percentile_vol_nutricula = calc_percentile(b_nutricula_angles(:,3),b_nutricula_angles(:,3));

%size versus mean angle and angle difference
figure
hold on
scatter(percentile_vol_pilosum,(s_pilosum_angles(:,1)+s_pilosum_angles(:,2))./2,"filled",'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560])
scatter(percentile_vol_minutum,(s_minutum_angles(:,1)+s_minutum_angles(:,2))./2,"filled",'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])
scatter(percentile_vol_nutricula,(b_nutricula_angles(:,1)+b_nutricula_angles(:,2))./2,"filled",'MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840])
xlabel('Chromosome Volume Percentile','FontSize', 24)
ylabel('$<\theta>$ [Degrees]','FontSize', 24, 'Interpreter','latex')
yline(90, '--', '90$^{\circ}$','FontSize', 24,'Interpreter','latex')
ylim([0 180])
xlim([0 100])
lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northeast', 'FontAngle', 'italic');
legend boxoff
ax = gca;
ax.FontSize = 16;

figure
hold on
scatter(percentile_vol_pilosum,s_pilosum_angles(:,1)-s_pilosum_angles(:,2),"filled",'MarkerEdgeColor',[0.4940 0.1840 0.5560],'MarkerFaceColor',[0.4940 0.1840 0.5560])
scatter(percentile_vol_minutum,s_minutum_angles(:,1)-s_minutum_angles(:,2),"filled",'MarkerEdgeColor',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250])
scatter(percentile_vol_nutricula,b_nutricula_angles(:,1)-b_nutricula_angles(:,2),"filled",'MarkerEdgeColor',[0.6350 0.0780 0.1840],'MarkerFaceColor',[0.6350 0.0780 0.1840])
xlabel('Chromosome Volume Percentile','FontSize', 24)
ylabel('$\Delta\theta$ [Degrees]','FontSize', 24, 'Interpreter','latex')
yline(0, '--', '0$^{\circ}$','FontSize', 24,'Interpreter','latex')
ylim([-90 90])
xlim([0 100])
lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northeast', 'FontAngle', 'italic');
legend boxoff
ax = gca;
ax.FontSize = 16;

function percentiles = calc_percentile(total_list,subset)
% Sort the total list to find the ranks
sorted_list = sort(total_list);

% Initialize array to store percentiles
percentiles = zeros(size(subset,1),1);

% Calculate percentile for each number in the subset
for i = 1:size(subset,1)
    % Find the rank of the current element in the sorted list
    rank = find(sorted_list == subset(i));
    
    % The percentile is the rank divided by the total number of elements
    percentiles(i) = (rank(1) - 1) / (size(sorted_list,1) - 1) * 100;
end
end


