%Lucas Philipp
clc;clear;

% Iterate through folder
ROI_start = 2;
ROI_end = 162;

for i = ROI_start:1:ROI_end
    I = tiffreadVolume(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/ROI %i.tiff",i));
    %I = tiffreadVolume(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium minutum/minutum cell 3 chromosomes/ROI %i.tiff",i));
    %Permute and project volume
    imwrite(project_back(I),sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Brandtodinium nutricula/Brandtodinium FL cell 1 chromosomes/back/projected_ROI_back%i.png",i));
end

% Function to permute and project
function [I3] = project_back(I)
%I2 = mean(I,1);
I2 = mean(I(round(size(I,1)/2):1:end,:,:),1);
I2 = permute(I2,[3 2 1]); %THIS CHANGES THE HANDEDNESS
I3 = (I2./max(I2)); %Normalize Image values
I3 = permute(I3,[2 1]); %THIS CHANGES THE HANDEDNESS BACK TO THE ORIGINAL HANDEDNESS
I3 = imrotate(I3,90);
end