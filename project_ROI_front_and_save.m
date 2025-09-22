%Lucas Philipp
%clc;clear;

% Iterate through folder
ROI_start = 1;
ROI_end = 130;

for i = ROI_start:1:ROI_end
    I = tiffreadVolume(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 3 Chromosomes 4nm voxels/ROI %i.tiff",i));
    %Permute and project volume
    imwrite(project_front(I),sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Ross Sea Dinoflagellate/Ross Sea Dinoflagellate Cell 3 Chromosomes 4nm voxels/front/projected_ROI_front%i.png",i));
end

% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder1/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder1.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder2/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder2.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder3/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder3.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder4/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder4.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder5/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder5.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder6/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder6.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder7/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder7.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder8/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder8.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder9/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder9.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder10/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder10.png');
% I = tiffreadVolume('/Users/lucasphilipp/Downloads/run3 se2 right cell edited/cylinder11/ROI 1.tiff');
% imwrite(project_front(I),'projected_ROI_front_cylinder11.png');

% Function to permute and project
function [I3] = project_front(I)
I2 = mean(I(1:1:round(size(I,1)/2),:,:),1);
I2 = permute(I2,[3 2 1]); %THIS CHANGES THE HANDEDNESS
I3 = (I2./max(I2)); %Normalize Image values
I3 = permute(I3,[2 1]); %THIS CHANGES THE HANDEDNESS BACK TO THE ORIGINAL HANDEDNESS
I3 = imrotate(I3,90);
end