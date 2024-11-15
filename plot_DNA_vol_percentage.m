
% %average genome volumes
% 
% microadriaticum cell 1 2.1169e+09
% 
% pilosum cell 1 2.5083e+09
% pilosum cell 2 2.7563e+09
% pilosum cell 3 2.9535e+09
% 
% minutum cell 1 2.7233e+09
% minutum cell 2 2.7608e+09
% minutum cell 3 2.3048e+09
% 
% brandotinium nutricula cell 1 1.5224e+10
% brandotinium nutricula cell 2 1.8176e+09
% brandotinium nutricula cell 3 1.6883e+10

% %average chromosome volumes
% 
% microadriaticum cell 1 2.4615e+07
% 
% pilosum cell 1 2.5336e+07
% pilosum cell 2 2.8126e+07
% pilosum cell 3 3.1420e+07
% 
% minutum cell 1 10.893e+07
% minutum cell 2  9.5201e+07
% minutum cell 3 7.4348e+07
% 
% brandotinium nutricula cell 1 1.2082e+08
% brandotinium nutricula cell 2 1.8738e+07
% brandotinium nutricula cell 3 5.9869e+07
% 
% %DNA vol/Chromatin vol

microadriaticum_DNA_Chromatin_vol_percentage = [55.7421; 0; 0];

pilosum_DNA_Chromatin_vol_percentage = [79.3363; 72.1970; 67.3773];

minutum_DNA_Chromatin_vol_percentage = [55.0797; 54.3317; 65.0815];

kawagutii_DNA_Chromatin_vol_percentage = [2.5; 2.1; 2.5];

%ADD COLOUR GRADIENT FOR DIFFERENT CELLS
%FaceColor = [0.4660 0.6740 0.1880]*1.45)
%FaceColor = [0.4660 0.6740 0.1880])
%FaceColor = [0.4660 0.6740 0.1880]/1.45)

% h=bar(x,'grouped') where x is an n by m matrix, will create a bar chart with n groups that each have m bars. h will be a 1xm matrix of bar handles.  
% To apply a color to each group, create a n by 3 matrix of RGB color values that defines the color for each group.  Then apply the color matrix to the CData property of each bar object (R2017b or later). You'll also need to set FaceColor to flat.  

data = [microadriaticum_DNA_Chromatin_vol_percentage'; pilosum_DNA_Chromatin_vol_percentage'; minutum_DNA_Chromatin_vol_percentage'; kawagutii_DNA_Chromatin_vol_percentage'];
figure
b = bar(data, 'grouped');

colors = [0 0.4470 0.7410; 0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880];
b(1).FaceColor = 'flat'; 
b(1).CData = colors;

b(2).FaceColor = 'flat'; 
b(2).CData = colors;

b(3).FaceColor = 'flat'; 
b(3).CData = colors;

set(gca,'XTickLabel',{'Symbiodinium microadriaticum','Symbiodinium pilosum','Symbiodinium minutum','Symbiodinium kawagutii'});
%lgd=legend({'Symbiodinium microadriaticum','Symbiodinium pilosum','Symbiodinium minutum','Symbiodinium kawagutii'},'Location', 'northwest', 'FontAngle', 'italic');
%legend boxoff
ylabel('(DNA Volume/Chromatin Volume)*100%','FontSize', 24)
ylim([0 100])
ax = gca;
ax.FontSize = 16;