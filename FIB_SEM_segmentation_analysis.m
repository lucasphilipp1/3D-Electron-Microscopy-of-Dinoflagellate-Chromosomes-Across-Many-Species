clc
clear
asphericity_tensor_cell1 = [];
ratio_P_axis_1_2_cell1 = [];
ratio_P_axis_1_3_cell1 = [];
ratio_P_axis_2_3_cell1 = [];

volume_cell1 = []; %in units number of voxels (16x16x16nm^3)

asphericity_tensor_cell2 = [];
ratio_P_axis_1_2_cell2 = [];
ratio_P_axis_1_3_cell2 = [];
ratio_P_axis_2_3_cell2 = [];

volume_cell2 = []; %in units number of voxels (16x16x16nm^3)

asphericity_tensor_cell3 = [];
ratio_P_axis_1_2_cell3 = [];
ratio_P_axis_1_3_cell3 = [];
ratio_P_axis_2_3_cell3 = [];

volume_cell3 = []; %in units number of voxels (16x16x16nm^3)

for j = 1:1:3
    if j == 1
    for i = 4:1:102

        V_format_PCA = [];

        V = tiffreadVolume(append(sprintf('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell %i',j), sprintf('/ROI %i.tiff',i)));

        %compute chromosome volume
        V=V./255;
        volume_cell1 = [volume_cell1; sum(sum(sum(V)))*16]; %in units number of voxels (16x16x16nm^3)

        stats = regionprops3(V,'PrincipalAxisLength','EigenValues');

        PC_cell1 = stats.PrincipalAxisLength;
        ratio_P_axis_1_2_cell1 = [ratio_P_axis_1_2_cell1; PC_cell1(1)/PC_cell1(2)];
        ratio_P_axis_1_3_cell1 = [ratio_P_axis_1_3_cell1; PC_cell1(1)/PC_cell1(3)];
        ratio_P_axis_2_3_cell1 = [ratio_P_axis_2_3_cell1; PC_cell1(2)/PC_cell1(3)];

        %calculate moment of inertia tensor
        %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
        eigenvalues = stats.EigenValues{1};

        asphericity_tensor_cell1 = [asphericity_tensor_cell1; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        i
    end
    elseif j ==2
        for i = 22:1:124

        V_format_PCA = [];

        V = tiffreadVolume(append(sprintf('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell %i',j), sprintf('/ROI %i.tiff',i)));

        %compute chromosome volume
        V=V./255;
        volume_cell2 = [volume_cell2; sum(sum(sum(V)))*16]; %in units number of voxels (16x16x16nm^3)

        stats = regionprops3(V,'PrincipalAxisLength','EigenValues');

        PC_cell2 = stats.PrincipalAxisLength;
        ratio_P_axis_1_2_cell2 = [ratio_P_axis_1_2_cell2; PC_cell2(1)/PC_cell2(2)];
        ratio_P_axis_1_3_cell2 = [ratio_P_axis_1_3_cell2; PC_cell2(1)/PC_cell2(3)];
        ratio_P_axis_2_3_cell2 = [ratio_P_axis_2_3_cell2; PC_cell2(2)/PC_cell2(3)];

        %calculate moment of inertia tensor
        %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
        eigenvalues = stats.EigenValues{1};

        asphericity_tensor_cell2 = [asphericity_tensor_cell2; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        i
        end
    elseif j==3
        for i = 27:1:124

        V_format_PCA = [];

        V = tiffreadVolume(append(sprintf('/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell %i',j), sprintf('/ROI %i.tiff',i)));

        %compute chromosome volume
        V=V./255;
        volume_cell3 = [volume_cell3; sum(sum(sum(V)))*16]; %in units number of voxels (16x16x16nm^3)

        stats = regionprops3(V,'PrincipalAxisLength','EigenValues');

        PC_cell3 = stats.PrincipalAxisLength;
        ratio_P_axis_1_2_cell3 = [ratio_P_axis_1_2_cell3; PC_cell3(1)/PC_cell3(2)];
        ratio_P_axis_1_3_cell3 = [ratio_P_axis_1_3_cell3; PC_cell3(1)/PC_cell3(3)];
        ratio_P_axis_2_3_cell3 = [ratio_P_axis_2_3_cell3; PC_cell3(2)/PC_cell3(3)];

        %calculate moment of inertia tensor
        %Following: Chu, X., & Wang, J. (2023). Quantifying the large-scale chromosome structural dynamics during the mitosis-to-G1 phase transition of cell cycle. Open Biology, 13(11), 230175.
        eigenvalues = stats.EigenValues{1};

        asphericity_tensor_cell3 = [asphericity_tensor_cell3; (3/(2*sum(eigenvalues)^2))*((eigenvalues(1)-mean(eigenvalues))^2+(eigenvalues(2)-mean(eigenvalues))^2+(eigenvalues(3)-mean(eigenvalues))^2)];

        i
    end
    end 
end

figure
hold on
[counts1, binCenters1] = hist(asphericity_tensor_cell1, 10);
[counts2, binCenters2] = hist(asphericity_tensor_cell2, 10);
[counts3, binCenters3] = hist(asphericity_tensor_cell3, 10);
plot(binCenters1, counts1, 'r-', 'LineWidth',2);
plot(binCenters2, counts1, 'g-', 'LineWidth',2);
plot(binCenters3, counts1, 'b-', 'LineWidth',2);
xlabel('Chromosome Asphericity','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;
xlim([0 0.45])
% Put up legend.
legend1 = sprintf('Cell 1');
legend2 = sprintf('Cell 2');
legend3 = sprintf('Cell 3');
lgd=legend({legend1, legend2, legend3});
legend boxoff
lgd.FontSize = 20;
title('Symbiodinium Clade A')

figure
hold on
[counts1, binCenters1] = hist(volume_cell1, 10);
[counts2, binCenters2] = hist(volume_cell2, 10);
[counts3, binCenters3] = hist(volume_cell3, 10);
plot(binCenters1, counts1, 'r-', 'LineWidth',2);
plot(binCenters2, counts2, 'g-', 'LineWidth',2);
plot(binCenters3, counts3, 'b-', 'LineWidth',2);
xlabel('Chromosome Volume [nm^3]','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;
% Put up legend.
legend1 = sprintf('Cell 1');
legend2 = sprintf('Cell 2');
legend3 = sprintf('Cell 3');
lgd=legend({legend1, legend2, legend3});
legend boxoff
lgd.FontSize = 20;
title('Symbiodinium Clade A')

figure
hold on
scatter(volume_cell1, asphericity_tensor_cell1,'filled','o','MarkerFaceColor','red','MarkerEdgeColor','red');
scatter(volume_cell2, asphericity_tensor_cell2,'filled','o','MarkerFaceColor','green','MarkerEdgeColor','green');
scatter(volume_cell3, asphericity_tensor_cell3,'filled','o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
xlabel('Chromosome Volume [nm^3]','FontSize', 24)
ylabel('Chromosome Asphericity','FontSize', 24)
ax = gca;
ax.FontSize = 16;
% Put up legend.
legend1 = sprintf('Cell 1');
legend2 = sprintf('Cell 2');
legend3 = sprintf('Cell 3');
lgd=legend({legend1, legend2, legend3},'Location','northwest');
legend boxoff
lgd.FontSize = 20;
title('Symbiodinium Clade A')

figure
hold on
[counts1, binCenters1] = hist(ratio_P_axis_1_2_cell1, 10);
[counts2, binCenters2] = hist(ratio_P_axis_1_3_cell1, 10);
[counts3, binCenters3] = hist(ratio_P_axis_2_3_cell1, 10);

[counts4, binCenters4] = hist(ratio_P_axis_1_2_cell2, 10);
[counts5, binCenters5] = hist(ratio_P_axis_1_3_cell2, 10);
[counts6, binCenters6] = hist(ratio_P_axis_2_3_cell2, 10);

[counts7, binCenters7] = hist(ratio_P_axis_1_2_cell3, 10);
[counts8, binCenters8] = hist(ratio_P_axis_1_3_cell3, 10);
[counts9, binCenters9] = hist(ratio_P_axis_2_3_cell3, 10);

plot(binCenters1, counts1, 'r-', 'LineWidth',2);
hold on;
plot(binCenters2, counts2, 'r-','LineWidth',2);
plot(binCenters3, counts3, 'r-','LineWidth',2);

plot(binCenters4, counts4, 'g-', 'LineWidth',2);
plot(binCenters5, counts5, 'g-','LineWidth',2);
plot(binCenters6, counts6, 'g-','LineWidth',2);

plot(binCenters7, counts7, 'b-', 'LineWidth',2);
plot(binCenters8, counts8, 'b-','LineWidth',2);
plot(binCenters9, counts9, 'b-','LineWidth',2);

xlabel('Asphericity (Ratio of Principle Components)','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;
xlim([1 3])
% Put up legend.
legend1 = sprintf('Cell 1');
legend2 = sprintf('Cell 2');
legend3 = sprintf('Cell 3');
lgd=legend({legend1, '', '', legend2, '', '', legend3, '', ''});
legend boxoff
lgd.FontSize = 20;
title('Symbiodinium Clade A')
