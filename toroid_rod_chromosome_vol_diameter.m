clc
clear

toroid_vol = [];
toroid_diameter = [];

rod_vol = [];
rod_diameter = [];

chromosome_vol_Brandtodinium_nutricula = [];
chromosome_vol_Symbiodinium_pilosum = [];
chromosome_vol_Symbiodinium_microadriaticum = [];
chromosome_vol_Symbiodinium_minutum = [];

cdBN = [];
cdSP = [];
cdSMic = [];
cdSMin = [];

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Brandtodinium nutricula/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        if ii == 1
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*10^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            cdBN = [cdBN; 8*s.PrincipalAxisLength(1) 8*s.PrincipalAxisLength(2) 8*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 2
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*4^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            cdBN = [cdBN; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 3
            chromosome_vol_Brandtodinium_nutricula = [chromosome_vol_Brandtodinium_nutricula; sum(I,"all")*2^3 ii]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            cdBN = [cdBN; 2*s.PrincipalAxisLength(1) 2*s.PrincipalAxisLength(2) 2*s.PrincipalAxisLength(3)]; %in nm
        end
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium microadriaticum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_microadriaticum = [chromosome_vol_Symbiodinium_microadriaticum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        cdSMic = [cdSMic; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium minutum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_minutum = [chromosome_vol_Symbiodinium_minutum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        cdSMin = [cdSMin; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Chromosomes/Symbiodinium pilosum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        chromosome_vol_Symbiodinium_pilosum = [chromosome_vol_Symbiodinium_pilosum; sum(I,"all")*4^3 ii]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        cdSP = [cdSP; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

%%% Toroids

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Toroids/Brandtodinium nutricula/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        if ii == 1
            toroid_vol = [toroid_vol; sum(I,"all")*10^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            toroid_diameter = [toroid_diameter; 8*s.PrincipalAxisLength(1) 8*s.PrincipalAxisLength(2) 8*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 2
            toroid_vol = [toroid_vol; sum(I,"all")*4^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            toroid_diameter = [toroid_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 3
            toroid_vol = [toroid_vol; sum(I,"all")*2^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            toroid_diameter = [toroid_diameter; 2*s.PrincipalAxisLength(1) 2*s.PrincipalAxisLength(2) 2*s.PrincipalAxisLength(3)]; %in nm
        end
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Toroids/Symbiodinium microadriaticum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        toroid_vol = [toroid_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        toroid_diameter = [toroid_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Toroids/Symbiodinium minutum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        toroid_vol = [toroid_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        toroid_diameter = [toroid_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Toroids/Symbiodinium pilosum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        toroid_vol = [toroid_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        toroid_diameter = [toroid_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

%%% Rods

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Rods/Brandtodinium nutricula/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        if ii == 1
            rod_vol = [rod_vol; sum(I,"all")*10^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            rod_diameter = [rod_diameter; 8*s.PrincipalAxisLength(1) 8*s.PrincipalAxisLength(2) 8*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 2
            rod_vol = [rod_vol; sum(I,"all")*4^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            rod_diameter = [rod_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
        elseif ii == 3
            rod_vol = [rod_vol; sum(I,"all")*2^3]; %in nm^3
            s = regionprops3(I,"PrincipalAxisLength");
            rod_diameter = [rod_diameter; 2*s.PrincipalAxisLength(1) 2*s.PrincipalAxisLength(2) 2*s.PrincipalAxisLength(3)]; %in nm
        end
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Rods/Symbiodinium microadriaticum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        rod_vol = [rod_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        rod_diameter = [rod_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Rods/Symbiodinium minutum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        rod_vol = [rod_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        rod_diameter = [rod_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

D = '/Users/lucasphilipp/Desktop/Research/Weber/Dinoflagellate FIB-SEM Data/Slice & View/Rods/Symbiodinium pilosum/';
S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
for ii = 1:numel(N)
    T = dir(fullfile(D,N{ii},'*.tiff')); % improve by specifying the file extension.
    C = {T(~[T.isdir]).name}; % files in subfolder.
    for jj = 1:numel(C)
        F = fullfile(D,N{ii},C{jj});
        I = tiffreadVolume(F);
        I=I./max(I);
        rod_vol = [rod_vol; sum(I,"all")*4^3]; %in nm^3
        s = regionprops3(I,"PrincipalAxisLength");
        rod_diameter = [rod_diameter; 4*s.PrincipalAxisLength(1) 4*s.PrincipalAxisLength(2) 4*s.PrincipalAxisLength(3)]; %in nm
    end
end

figure
hold on
histogram(toroid_diameter(:,3), FaceColor = [0.25 0.25 0.25])
histogram((toroid_diameter(:,1)+toroid_diameter(:,2))/2, FaceColor = [0.75 0.75 0.75])
ax = gca;
ax.FontSize = 16;
lgd=legend({'Toroid Thickness','Toroid Diameter'});
legend boxoff
xlabel('DNA Toroid Dimensions [nm]')
xlim([0 1000])
ylabel('Count')

figure
hold on
histogram((rod_diameter(:,2)+rod_diameter(:,3))/2, FaceColor = [0.25 0.25 0.25])
histogram(rod_diameter(:,1), FaceColor = [0.75 0.75 0.75])
ax = gca;
ax.FontSize = 16;
lgd=legend({'Rod Diameter','Rod Length'});
legend boxoff
xlabel('DNA Rod Dimensions [nm]')
xlim([0 1000])
ylabel('Count')

figure
hold on
histogram(rod_vol, FaceColor = [0 0 0], FaceAlpha = 0.8)
histogram(toroid_vol, FaceColor = [1 1 1])
ax = gca;
ax.FontSize = 16;
lgd=legend({'Rods','Toroids'});
legend boxoff
xlabel('Volume [nm^3]')
ylabel('Count')

edges=linspace(0,10^8,20);
figure
hold on
histogram(chromosome_vol_Symbiodinium_microadriaticum(find(chromosome_vol_Symbiodinium_microadriaticum(:,2)==1),1), edges, FaceColor = [0 0.4470 0.7410])
ax = gca;
ax.FontSize = 24;
lgd=legend({'Cell 1'});
%lgd=legend({'Cell 1', 'Cell 2', 'Cell 3'});
%xlim([0 6*10^8])
xlim([0 10*10^7])
legend boxoff
xlabel('Chromosome Volume [nm^3]')
ylabel('Count')
title('Symbiodinium microadriaticum')

edges=linspace(0,10^8,20);
figure
hold on
histogram(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==1),1), edges, FaceColor = [0.4940 0.1840 0.5560])
histogram(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==2),1), edges, FaceColor = [0.4940 0.1840 0.5560]./1.5)
histogram(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==3),1), edges, FaceColor = [0.4940 0.1840 0.5560]./2)
ax = gca;
ax.FontSize = 24;
lgd=legend({'Cell 1', 'Cell 2', 'Cell 3'});
%xlim([0 6*10^8])
xlim([0 10*10^7])
legend boxoff
xlabel('Chromosome Volume [nm^3]')
ylabel('Count')
title('Symbiodinium pilosum')

edges=linspace(0,40*10^7,20);
figure
hold on
histogram(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==1),1), edges, FaceColor = [0.9290 0.6940 0.1250])
histogram(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==2),1), edges, FaceColor = [0.9290 0.6940 0.1250]./1.5)
histogram(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==3),1), edges, FaceColor = [0.9290 0.6940 0.1250]./2)
ax = gca;
ax.FontSize = 24;
lgd=legend({'Cell 1', 'Cell 2', 'Cell 3'});
xlim([0 6*10^8])
legend boxoff
xlabel('Chromosome Volume [nm^3]')
ylabel('Count')
title('Symbiodinium minutum')

edges=linspace(0,6*10^8,20);
figure
hold on
histogram(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==1),1), edges, FaceColor = [0.6350 0.0780 0.1840])
histogram(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==2),1), edges, FaceColor = [0.6350 0.0780 0.1840]./1.5)
histogram(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==3),1), edges, FaceColor = [0.6350 0.0780 0.1840]./2)
ax = gca;
ax.FontSize = 24;
lgd=legend({'Cell 1', 'Cell 2', 'Cell 3'});
xlim([0 6*10^8])
legend boxoff
xlabel('Chromosome Volume [nm^3]')
ylabel('Count')
title('Brandtodinium nutricula')

edges=linspace(1,5,100);
figure
hold on
histogram(((cdSMic(:,1)./cdSMic(:,2))+(cdSMic(:,1)./cdSMic(:,3))+(cdSMic(:,2)./cdSMic(:,3)))./3, edges, FaceColor = [0 0.4470 0.7410])
histogram(((cdSP(:,1)./cdSP(:,2))+(cdSP(:,1)./cdSP(:,3))+(cdSP(:,2)./cdSP(:,3)))./3, edges, FaceColor = [0.4940 0.1840 0.5560])
histogram(((cdSMin(:,1)./cdSMin(:,2))+(cdSMin(:,1)./cdSMin(:,3))+(cdSMin(:,2)./cdSMin(:,3)))./3, edges, FaceColor = [0.9290 0.6940 0.1250])
histogram(((cdBN(:,1)./cdBN(:,2))+(cdBN(:,1)./cdBN(:,3))+(cdBN(:,2)./cdBN(:,3)))./3, edges, FaceColor = [0.6350 0.0780 0.1840])
ax = gca;
ax.FontSize = 24;
lgd=legend({'Symbiodinium microadriaticum','Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'});
legend boxoff
xlabel('Chromosome Asphericity')
ylabel('Count')

%genome volumes
sum(chromosome_vol_Symbiodinium_microadriaticum(find(chromosome_vol_Symbiodinium_microadriaticum(:,2)==1),1))

sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==1),1))
sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==2),1))
sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==3),1))

sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==1),1))
sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==2),1))
sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==3),1))

sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==1),1))
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==2),1))
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==3),1))

%average chromosome volumes
sum(chromosome_vol_Symbiodinium_microadriaticum(find(chromosome_vol_Symbiodinium_microadriaticum(:,2)==1),1))/size(find(chromosome_vol_Symbiodinium_microadriaticum(:,2)==1),1)

sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==1),1))/size(find(chromosome_vol_Symbiodinium_pilosum(:,2)==1),1)
sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==2),1))/size(find(chromosome_vol_Symbiodinium_pilosum(:,2)==2),1)
sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==3),1))/size(find(chromosome_vol_Symbiodinium_pilosum(:,2)==3),1)

sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==1),1))/size(find(chromosome_vol_Symbiodinium_minutum(:,2)==1),1)
sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==2),1))/size(find(chromosome_vol_Symbiodinium_minutum(:,2)==2),1)
sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==3),1))/size(find(chromosome_vol_Symbiodinium_minutum(:,2)==3),1)

sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==1),1))/size(find(chromosome_vol_Brandtodinium_nutricula(:,2)==1),1)
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==2),1))/size(find(chromosome_vol_Brandtodinium_nutricula(:,2)==2),1)
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==3),1))/size(find(chromosome_vol_Brandtodinium_nutricula(:,2)==3),1)

%DNA vol/genome vol*100%
%genome assemblies estimate: microadriaticum 1.18 Gbp, pilosum 1.99 Gbp, minutum 1.50 Gbp, kawagutii 1.10 Gbp
(1.18*10^9)*100/sum(chromosome_vol_Symbiodinium_microadriaticum(find(chromosome_vol_Symbiodinium_microadriaticum(:,2)==1),1))

(1.99*10^9)*100/sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==1),1))
(1.99*10^9)*100/sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==2),1))
(1.99*10^9)*100/sum(chromosome_vol_Symbiodinium_pilosum(find(chromosome_vol_Symbiodinium_pilosum(:,2)==3),1))

(1.50*10^9)*100/sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==1),1))
(1.50*10^9)*100/sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==2),1))
(1.50*10^9)*100/sum(chromosome_vol_Symbiodinium_minutum(find(chromosome_vol_Symbiodinium_minutum(:,2)==3),1))

sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==1),1))
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==2),1))
sum(chromosome_vol_Brandtodinium_nutricula(find(chromosome_vol_Brandtodinium_nutricula(:,2)==3),1))
