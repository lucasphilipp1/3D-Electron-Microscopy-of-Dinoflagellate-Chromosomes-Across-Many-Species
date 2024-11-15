clc
clear
for i = 2:1:35
    V = tiffreadVolume(sprintf('ROI %i.tiff',i));
    % create viewer
    viewer = viewer3d(BackgroundColor="white",Lighting="on",BackgroundGradient="off");
    % volshow
    volshow(V,Parent=viewer);
    pause(5) % pause to give viewer time to fully open
    % capture frame
    frame = getframe(viewer.Parent);
    imwrite(frame2im(frame),append(sprintf('ROIvolshow%i',i),'.png'))
    %viewer.clear % clear previous   
end
