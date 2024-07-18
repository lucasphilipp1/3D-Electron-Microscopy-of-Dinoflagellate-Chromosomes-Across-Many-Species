s_pilosum_angles = [];
b_nutricula_angles = [];
s_minutum_angles =[];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 1/Clade_A_cell_1_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 2/Clade_A_cell_2_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium pilosum/Cell 3/Clade_A_cell_3_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
s_pilosum_angles = [s_pilosum_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/brandtodinium_FL_cell_1_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium cell 2 stack chromosomes/brandtodinium_cell_stack_2_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 3 chromosomes/brandtodinium_FL_cell_3_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
b_nutricula_angles = [b_nutricula_angles; ImageJ_angles(:,7)];

ImageJ_angles = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum_chromosome_angles.csv");
ImageJ_angles(1,:)=[]; %remove header
s_minutum_angles = [s_minutum_angles; ImageJ_angles(:,7)];

figure
hold on
histogram(s_pilosum_angles,15, 'EdgeColor', [0.4940 0.1840 0.5560], 'FaceColor', [0.4940 0.1840 0.5560])
histogram(s_minutum_angles,15, 'EdgeColor', [0.9290 0.6940 0.1250], 'FaceColor', [0.9290 0.6940 0.1250])
histogram(b_nutricula_angles,15, 'EdgeColor', [0.6350 0.0780 0.1840], 'FaceColor', [0.6350 0.0780 0.1840])
xlim([0 180])
xlabel('Angle [Degrees]','FontSize', 24)
ylabel('Count','FontSize', 24)
lgd=legend({'Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'},'Location', 'northwest', 'FontAngle', 'italic');
legend boxoff
ax = gca;
ax.FontSize = 16;


