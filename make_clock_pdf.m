clc
clear

ROI_start = 2;
ROI_end = 115;

Images_front = {ROI_end,1};
Images_back = {ROI_end,1};

Images_front_clock = {ROI_end,1};
Images_back_clock = {ROI_end,1};

ImageJ_angles = readtable("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ET cell 1.csv");
clock_info = [str2double(extract(ImageJ_angles.Label, digitsPattern)) ImageJ_angles.Angle]; %back first, front second

for i = ROI_start:1:ROI_end
    Images_front{i} = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes/front/projected_ROI_front%i.png",i));
    Images_back{i} = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes/back/projected_ROI_back%i.png",i));
end

for j = 1:1:size(clock_info,1)
    i=clock_info(j,1);
    Images_front_clock{i} = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes/check/angle_check_front_ROI_%i.jpg",i));
    Images_back_clock{i} = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ensiculifera tyrrhenica/ensiculifera tyrrhenica chromosomes/check/angle_check_back_ROI_%i.jpg",i));
end

for i = ROI_start:1:ROI_end
    if size(find(clock_info(:,1)==i),1)>0 %angle extracted

        angles = clock_info(find(clock_info(:,1)==i),2);

        figure;
        tiledlayout(1,4)
        nexttile
        imshow(Images_front{i});
        title('Front Half');
        nexttile
        imshow(Images_front_clock{i});
        title('Extracted Angle');
        hold off
        sgtitle(sprintf('ROI %i',i));
        nexttile
        imshow(Images_back{i});
        title('Back Half');
        nexttile
        imshow(Images_back_clock{i});
        title('Extracted Angle');
        hold off
        sgtitle(sprintf('ROI %i',i));
    else %no angle extracted
        figure;
        tiledlayout(1,2)
        nexttile
        imshow(Images_front{i});
        title('Front Half');
        nexttile
        imshow(Images_back{i});
        title('Back Half');
        sgtitle(sprintf('ROI %i',i));
    end
    i
end

% append each of the figures to output.pdf
for i=1:ROI_end-ROI_start+1
    exportgraphics(figure(i), 'Ensiculifera tyrrhenica Cell 1 chromosome angles.pdf', 'Append', true);
end

% close all
% 
% orientation_front = [];
% orientation_back = [];
% 
% % correlate image with itself displaced in a direction
% for i = ROI_start:1:ROI_end
%     output = zeros(size(Images_front{i}));
%     output=stics(Images_front{i},1);
%     mask = output>(max(output,[],'all')-min(output,[],'all'))/2;
%     rp = regionprops(mask,'Orientation');
%     orientation_front = [orientation_front; rp.Orientation];
% 
%     output = zeros(size(Images_back{i}));
%     output=stics(Images_back{i},1);
%     mask = output>(max(output,[],'all')-min(output,[],'all'))/2;
%     rp = regionprops(mask,'Orientation');
%     orientation_back = [orientation_back; rp.Orientation];
% end