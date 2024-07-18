clear all;
close all;

min_branch_length = 3; %length in pixels
%keep track of downsampling and of voxel size
voxel_edge_dimension = 4; %in nm
downsampling_factor = 4;

% load example binary skeleton image
vol = tiffreadVolume('/Users/lucasphilipp/Downloads/kawagutii run2 middle cell downsample by 4 ROI.tiff');
VolLogical = imbinarize(vol);
skel = Skeleton3D(VolLogical);

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node,link] = Skel2Graph3D(skel,min_branch_length);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)

    wl = wl_new;   
    
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

     wl_new = sum(cellfun('length',{node2.links}));

end;

% display result
figure();
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end;
        if(node2(i).ep==1)
            col='k';
        end;

        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);

scale_factor=voxel_edge_dimension*downsampling_factor;
A2_full_cell1 = full(A2);
A2_full_cell1 = scale_factor.*A2_full_cell1;

struct_node2 = struct2cell(node2);
connectivity_cell1 = cellfun(@numel, struct_node2(2,:,:));

% load example binary skeleton image
vol = tiffreadVolume('/Users/lucasphilipp/Downloads/kawagutii run2 top right cell downsample by 4 ROI.tiff');
VolLogical = imbinarize(vol);
skel = Skeleton3D(VolLogical);

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node,link] = Skel2Graph3D(skel,min_branch_length);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)

    wl = wl_new;   
    
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

     wl_new = sum(cellfun('length',{node2.links}));

end;

% display result
figure();
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end;
        if(node2(i).ep==1)
            col='k';
        end;

        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);

scale_factor=voxel_edge_dimension*downsampling_factor;
A2_full_cell2 = full(A2);
A2_full_cell2 = scale_factor.*A2_full_cell2;

struct_node2 = struct2cell(node2);
connectivity_cell2 = cellfun(@numel, struct_node2(2,:,:));

% load example binary skeleton image
vol = tiffreadVolume('/Users/lucasphilipp/Downloads/kawagutii run3 right cell downsample by 4 ROI.tiff');
VolLogical = imbinarize(vol);
skel = Skeleton3D(VolLogical);

w = size(skel,1);
l = size(skel,2);
h = size(skel,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node,link] = Skel2Graph3D(skel,min_branch_length);

% total length of network
wl = sum(cellfun('length',{node.links}));

skel2 = Graph2Skel3D(node,link,w,l,h);
[~,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

% calculate new total length of network
wl_new = sum(cellfun('length',{node2.links}));

% iterate the same steps until network length changed by less than 0.5%
while(wl_new~=wl)

    wl = wl_new;   
    
     skel2 = Graph2Skel3D(node2,link2,w,l,h);
     [A2,node2,link2] = Skel2Graph3D(skel2,min_branch_length);

     wl_new = sum(cellfun('length',{node2.links}));

end;

% display result
figure();
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end;
        if(node2(i).ep==1)
            col='k';
        end;

        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);

scale_factor=voxel_edge_dimension*downsampling_factor;
A2_full_cell3 = full(A2);
A2_full_cell3 = scale_factor.*A2_full_cell3;

struct_node2 = struct2cell(node2);
connectivity_cell3 = cellfun(@numel, struct_node2(2,:,:));

figure
hold on
histogram(A2_full_cell1(:),'BinLimits',[1,max(A2_full_cell1(:))])
histogram(A2_full_cell2(:),'BinLimits',[1,max(A2_full_cell2(:))])
histogram(A2_full_cell3(:),'BinLimits',[1,max(A2_full_cell3(:))])
xlabel('Edge Contour Lengths [nm]')
ylabel('Count')

figure
hold on
histogram(connectivity_cell1)
histogram(connectivity_cell2)
histogram(connectivity_cell3)
xlabel('Node Connectivity')
ylabel('Count')