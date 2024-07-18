clear all;
close all;

vol = tiffreadVolume('/Users/lucasphilipp/Downloads/kawagutii run3 right cell downsample by 4 ROI.tiff');

VolLogical = imbinarize(vol);

skel = Skeleton3D(VolLogical);

figure();
col=[.7 .7 .8];
hiso = patch(isosurface(VolLogical,0),'FaceColor',col,'EdgeColor','none');
hiso2 = patch(isocaps(VolLogical,0),'FaceColor',col,'EdgeColor','none');
axis equal;axis off;
lighting phong;
isonormals(VolLogical,hiso);
alpha(0.5);
set(gca,'DataAspectRatio',[1 1 1])
camlight;
hold on;
w=size(skel,1);
l=size(skel,2);
h=size(skel,3);
[x,y,z]=ind2sub([w,l,h],find(skel(:)));
plot3(y,x,z,'square','Markersize',4,'MarkerFaceColor','r','Color','r');            
set(gcf,'Color','white');
view(140,80)

volshow(skel)