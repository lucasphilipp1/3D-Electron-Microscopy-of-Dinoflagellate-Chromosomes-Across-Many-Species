%Eran Ittah and Lucas Philipp
clc;clear;

% Iterate through folder and Save all FFTs and residual Images
ROI_start = 22;
ROI_end = 124;

radial_domain_all_chromosomes = cell(ROI_end-ROI_start+1,1);
radialProfile_all_chromosomes = cell(ROI_end-ROI_start+1,1);

for i = ROI_start:1:ROI_end

    I = tiffreadVolume(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell 2/ROI %i.tiff",i));
    %Permute and project volume
    I2 = Permute_Project(I);
    I2_smoothed = smoothdata(I2,"gaussian","SmoothingFactor",0.5);
    Ires = abs(I2_smoothed-I2)./max(I2_smoothed-I2);

    Icrop = imcrop(Ires,[0.2*size(Ires,2) 0.2*size(Ires,1) 0.6*size(Ires,2) 0.6*size(Ires,1)]);
    Icrop(Icrop>prctile(Ires(:),60))=prctile(Ires(:),60); %truncate pixel intensity at percentile, prevents bright boundary
    imwrite(Icrop,sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell 2/input/residual_ROI_%i.png",i));

    Image_for_FT = Icrop;

    % Perform 2D FFTs
    fftOriginal = fft2(double(Image_for_FT));
    % Move center from (1,1) to (129, 129) (the middle of the matrix).
    shiftedFFT = fftshift(fftOriginal);
    scaledFFTr = 255 * mat2gray(real(shiftedFFT));
    shiftedFFTMagnitude = abs(shiftedFFT);
    FFT_final = log(shiftedFFTMagnitude);
 
    % Get the average radial profile
    [rows, columns, numberOfColorChannels] = size(Image_for_FT);
    midRow = rows/2+1;
    midCol = columns/2+1;
    maxRadius = ceil(sqrt((floor(min(size(Image_for_FT))/2))^2 + (floor(min(size(Image_for_FT))/2))^2));
    radialProfile = zeros(maxRadius, 1);
    count = zeros(maxRadius, 1);
    for col = 1 : columns
        for row = 1 : rows
            radius = sqrt((row - midRow) ^ 2 + (col - midCol) ^ 2);
            thisIndex = ceil(radius) + 1;
            if(radius<maxRadius-1)
                radialProfile(thisIndex) = radialProfile(thisIndex) + log(abs(shiftedFFTMagnitude(row, col))); %changed last term from shiftedFFTMagnitude(row, col)
                count(thisIndex) = count(thisIndex) + 1;
            end
        end
    end
    % Get average
    radialProfile = radialProfile ./ count;
    radialProfile_all_chromosomes{i-ROI_start+1} = radialProfile;
    radial_domain_all_chromosomes{i-ROI_start+1} = 1:1:size(radialProfile,1);
end

figure
hold on
for i = 1:1:ROI_end-ROI_start+1
    plot(radial_domain_all_chromosomes{i}, radialProfile_all_chromosomes{i},'Color', [0 0.4470 0.741 0.5])
end
ylabel("Radial Average of logFT Intensity")
xlabel("wavenumber [1/16nm]")
ax = gca;
ax.FontSize = 16;

% Function to permute and project
function [I3] = Permute_Project(I)
I2 = mean(I,1);
I2 = permute(I2,[3 2 1]);
%Normalize Matrix values
I3 = (I2./max(I2));
end