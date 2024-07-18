clc
clear

ROI_start = 22;
ROI_end = 124;

angles_for_histogram = [];
for i = ROI_start:1:ROI_end
    orientationJ_distribution = readmatrix("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell 2/output/"+sprintf("orientation_distribution_%i.csv",i-ROI_start));
    angle = orientationJ_distribution(find(orientationJ_distribution(:,2)==max(orientationJ_distribution(:,2))),1);
    angles_for_histogram = [angles_for_histogram; angle];

    Icheck = imread("/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell 2/input/"+sprintf("residual_ROI_%i.png",i));

    figure
    imshow(Icheck,[])
    hold
    %draw a line over image to check ridge angle measurement is correct
    %L is the length of the line
    %alpha is the angle of the line
    alpha = 270-angle; %orientationJ: +90 is up, -90 is down, 0 is right %code below: 0 down, 90 left, 180 up, 270 right
    L = size(Icheck,1)/cosd(alpha);
    x2=(L*cosd(90-alpha));
    y2=(L*sind(90-alpha));
    x3=(L*cosd(-(90-alpha)));
    y3=size(Icheck,1)+(L*sind(-(90-alpha)));
    %plot([-tand(alpha)*size(Image_for_FT,1)/2+size(Image_for_FT,2)/2 x2-tand(alpha)*size(Image_for_FT,1)/2+size(Image_for_FT,2)/2],[0 y2],'k','LineWidth',1)
    %plot([-tand(alpha)*size(Image_for_FT,1)/2+size(Image_for_FT,2)/2 x3-tand(alpha)*size(Image_for_FT,1)/2+size(Image_for_FT,2)/2],[size(Image_for_FT,1) y3],'k','LineWidth',1)
    stox=linspace((-tand(alpha)*size(Icheck,1)/2+size(Icheck,2)/2),(x3-tand(alpha)*size(Icheck,1)/2+size(Icheck,2)/2),100);
    stoy=linspace(size(Icheck,1),y3,100);
    if 90<alpha && alpha<270
        %text(stox(70),stoy(70), num2str(alpha),'Color','w','FontSize',2)
        plot([stox(50) stox(100)],[stoy(50) stoy(100)],'w','LineWidth',1)
    else
        %text(stox(30),stoy(30), num2str(alpha),'Color','w','FontSize',2)
        plot([stox(1) stox(50)],[stoy(1) stoy(50)],'w','LineWidth',1)
    end
    ax = gca;
    exportgraphics(ax,"/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Symbiodinium Clade A/Cell 2/check/"+sprintf('angle_check_ROI_%i.jpg',i),"Resolution",1000)
    pause(1)
    close
end

figure
histogram(angles_for_histogram)
xlabel('Angle [0 Degrees Means Perpendicular to Chromosome Long Axis]','FontSize', 24)
ylabel('Count','FontSize', 24)
ax = gca;
ax.FontSize = 16;

