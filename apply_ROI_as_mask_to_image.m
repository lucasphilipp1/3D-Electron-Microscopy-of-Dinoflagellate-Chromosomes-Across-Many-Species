clc
clear

isInverted = 0;

Image = double(tiffreadVolume('/Users/lucasphilipp/Downloads/Image.tiff'));
ROI = double(tiffreadVolume('/Users/lucasphilipp/Downloads/ROI.tiff'));

if isInverted == 0
    Mask = 1-(ROI./255);
    else
    Mask = ROI./255;
end

Masked_Image = Image.*Mask;

%Masked_Image = Image;
% Write flag: run once with iwrite = false if you want to see what you are writing first
iwrite = true;
% File name:
fname = 'Test_Masked_Image.tif';
for n = 1:size(Masked_Image, 3)
    if n == 1
        % First slice:
        imwrite(Masked_Image(:, :, n),fname)
    else
        % Subsequent slices:
        imwrite(Masked_Image(:, :, n),fname,'WriteMode','append');
    end
end

