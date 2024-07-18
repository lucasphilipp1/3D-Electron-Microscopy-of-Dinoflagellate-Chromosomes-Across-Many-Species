for i = 4:1:102
    V = tiffreadVolume(sprintf('/Users/lucasphilipp/Downloads/Cell1_ROI%i.tiff',i));

    output=max(V,[],[1 2]);
    output=reshape(output,[],1);

    zmin=min(find(output==max(output)));
    zmax=max(find(output==max(output)));

    output=max(V,[],[1 3]);
    output=reshape(output,[],1);

    ymin=min(find(output==max(output)));
    ymax=max(find(output==max(output)));

    output=max(V,[],[2 3]);
    output=reshape(output,[],1);

    xmin=min(find(output==max(output)));
    xmax=max(find(output==max(output)));

    V=V(xmin-10:xmax+10,ymin-10:ymax+10,zmin-10:zmax+10);
    volshow(V(xmin-10:xmax+10,ymin-10:ymax+10,zmin-10:zmax+10));

    for n = 1:size(V, 3)
        if n == 1
            % First slice:
            imwrite(V(:, :, n),sprintf('Cell1_ROI_cropped%i.tiff',i))
        else
            % Subsequent slices:
            imwrite(V(:, :, n),sprintf('Cell1_ROI_cropped%i.tiff',i),'WriteMode','append');
        end
    end
    i
end




