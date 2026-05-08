clc
clear

ROI_start = 2;
ROI_end = 187;

% %colours
% 'Symbiodinium microadriaticum'
% [0 0.4470 0.7410]
% 'Symbiodinium pilosum'
% [0.4940 0.1840 0.5560]
% 'Symbiodinium minutum'
% [0.9290 0.6940 0.1250]
% 'Brandtodinium nutricula'
% [0.6350 0.0780 0.1840]
% 'Ensiculifera tyrrhenica'
% [0.7 0.7 0.7]
% 'Ross Sea Dinoflagellate'
% [0 0 0]
% Symbiodinium kawagutii
% [0.4660 0.6740 0.1880]
% Crypthecodinium cohnii
% [1 0 1]

ImageJ_angles = readtable("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium Cohnii/CC cell 2.csv");
clock_info = [str2double(extract(ImageJ_angles.Label, digitsPattern)) ImageJ_angles.Angle]; %back first, front second

for i = ROI_start:1:ROI_end
    if size(find(clock_info(:,1)==i),1)>0 %angle extracted
        angles = clock_info(find(clock_info(:,1)==i),2);
        I_front = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/front/projected_ROI_front%i.png",i));
        figure
        imshow(I_front,[])
        hold
        %draw a line over image to check ridge angle measurement is correct
        %L is the length of the line
        %alpha is the angle of the line
        alpha = angles(2);
        L = size(I_front,1)/cosd(alpha);
        x2=(L*cosd(90-alpha));
        y2=(L*sind(90-alpha));
        x3=(L*cosd(-(90-alpha)));
        y3=size(I_front,1)+(L*sind(-(90-alpha)));
        stox=linspace((-tand(alpha)*size(I_front,1)/2+size(I_front,2)/2),(x3-tand(alpha)*size(I_front,1)/2+size(I_front,2)/2),5000);
        stoy=linspace(size(I_front,1),y3,5000);
        if alpha < 90
            % %colours
            % % 'Symbiodinium microadriaticum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0 0.4470 0.7410],'LineWidth',2)
            % % 'Symbiodinium pilosum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.4940 0.1840 0.5560],'LineWidth',2)
            % 'Symbiodinium minutum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.9290 0.6940 0.1250],'LineWidth',2)
            % 'Crypthecodinium cohnii'
            plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [1 0 1],'LineWidth',2)
            % % 'Brandtodinium nutricula'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.6350 0.0780 0.1840],'LineWidth',2)
            % 'Ensiculifera tyrrhenica'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.7 0.7 0.7],'LineWidth',2)
            % % 'Ross Sea Dinoflagellate'
            %plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.4 0.4 0.4],'LineWidth',2)
        else
            % %colours
            % % 'Symbiodinium microadriaticum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0 0.4470 0.7410],'LineWidth',2)
            % % 'Symbiodinium pilosum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.4940 0.1840 0.5560],'LineWidth',2)
            % 'Symbiodinium minutum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
            % 'Crypthecodinium cohnii'
            plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [1 0 1],'LineWidth',2)
            % % 'Brandtodinium nutricula'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.6350 0.0780 0.1840],'LineWidth',2)
            % 'Ensiculifera tyrrhenica'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.7 0.7 0.7],'LineWidth',2)
            % % 'Ross Sea Dinoflagellate'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.4 0.4 0.4],'LineWidth',2)
        end
        ax = gca;
        exportgraphics(ax,"/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/check/"+sprintf('angle_check_front_ROI_%i.jpg',i),"Resolution",1000)
        close

        I_back = imread(sprintf("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/back/projected_ROI_back%i.png",i));
        figure
        imshow(I_back,[])
        hold
        %draw a line over image to check ridge angle measurement is correct
        %L is the length of the line
        %alpha is the angle of the line
        alpha = angles(1);
        L = size(I_back,1)/cosd(alpha);
        x2=(L*cosd(90-alpha));
        y2=(L*sind(90-alpha));
        x3=(L*cosd(-(90-alpha)));
        y3=size(I_back,1)+(L*sind(-(90-alpha)));
        stox=linspace((-tand(alpha)*size(I_back,1)/2+size(I_back,2)/2),(x3-tand(alpha)*size(I_back,1)/2+size(I_back,2)/2),5000);
        stoy=linspace(size(I_back,1),y3,5000);
        if alpha < 90
            % %colours
            % % 'Symbiodinium microadriaticum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0 0.4470 0.7410],'LineWidth',2)
            % % 'Symbiodinium pilosum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.4940 0.1840 0.5560],'LineWidth',2)
            % 'Symbiodinium minutum'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.9290 0.6940 0.1250],'LineWidth',2)
            % 'Crypthecodinium cohnii'
            plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [1 0 1],'LineWidth',2)
            % % 'Brandtodinium nutricula'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.6350 0.0780 0.1840],'LineWidth',2)
            % % 'Ensiculifera tyrrhenica'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.7 0.7 0.7],'LineWidth',2)
            % % 'Ross Sea Dinoflagellate'
            % plot([stox(1) stox(2500)],[stoy(1) stoy(2500)],'Color', [0.4 0.4 0.4],'LineWidth',2)
        else
            % %colours
            % % 'Symbiodinium microadriaticum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0 0.4470 0.7410],'LineWidth',2)
            % % 'Symbiodinium pilosum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.4940 0.1840 0.5560],'LineWidth',2)
            % 'Symbiodinium minutum'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
            % 'Crypthecodinium cohnii'
            plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color',[1 0 1],'LineWidth',2)
            % % 'Brandtodinium nutricula'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.6350 0.0780 0.1840],'LineWidth',2)
            % % 'Ensiculifera tyrrhenica'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.7 0.7 0.7],'LineWidth',2)
            % % 'Ross Sea Dinoflagellate'
            % plot([stox(2500) stox(5000)],[stoy(2500) stoy(5000)],'Color', [0.4 0.4 0.4],'LineWidth',2)
        end
        ax = gca;
        exportgraphics(ax,"/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Crypthecodinium cohnii/Cell 2/check/"+sprintf('angle_check_back_ROI_%i.jpg',i),"Resolution",1000)
        close
    end
end
