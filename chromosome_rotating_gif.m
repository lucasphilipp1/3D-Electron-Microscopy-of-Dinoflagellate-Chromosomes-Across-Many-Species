V = tiffreadVolume('ROI 33.tiff');
viewer = viewer3d(BackgroundColor="white",Lighting="on",BackgroundGradient="off");
h = volshow(V,Parent=viewer);
viewer = h.Parent;
hFig = viewer.Parent;
drawnow
filename = "rotating_chromosome.gif";
sz = size(V);
center = sz/2 + 0.5;
viewer.CameraTarget = center;
numberOfFrames = 360;
vec = linspace(0,2*pi,numberOfFrames)';
dist = sqrt(sz(1)^2 + sz(2)^2 + sz(3)^2);
myPosition = center + ([cos(vec) sin(vec) zeros(size(vec))]*dist);
for idx = 1:length(vec)
    % Update the current view
    viewer.CameraPosition = myPosition(idx,:);
    % Capture the image using the getframe function
    I = getframe(hFig);
    [indI,cm] = rgb2ind(I.cdata,256);
    % Write the frame to the GIF file
    if idx==1
        % Do nothing. The first frame displays only the viewer, not the
        % volume.
    elseif idx == 2
        imwrite(indI,cm,filename,"gif",Loopcount=inf,DelayTime=0)
    else
        imwrite(indI,cm,filename,"gif",WriteMode="append",DelayTime=0)
    end
end