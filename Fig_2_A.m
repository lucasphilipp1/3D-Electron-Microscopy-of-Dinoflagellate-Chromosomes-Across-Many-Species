%microadriaticum, pilosum, minutum, brantodinium (3 cells each)
colors = [0 0.4470 0.7410; 0 0.4470 0.7410; 0 0.4470 0.7410; 0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; 0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; 0.9290 0.6940 0.1250; 0.6350 0.0780 0.1840; 0.6350 0.0780 0.1840; 0.6350 0.0780 0.1840];

chromosome_vol_Symbiodinium_microadriaticum

chromosome_vol_Symbiodinium_pilosum(:,2) = chromosome_vol_Symbiodinium_pilosum(:,2) + 1;

chromosome_vol_Symbiodinium_minutum(:,2) = chromosome_vol_Symbiodinium_minutum(:,2) + 4;

chromosome_vol_Brandtodinium_nutricula(:,2) = chromosome_vol_Brandtodinium_nutricula(:,2) + 7;

figure
hold on
histogram2(chromosome_vol_Symbiodinium_microadriaticum(:,1),chromosome_vol_Symbiodinium_microadriaticum(:,2),'FaceColor',[0 0.4470 0.7410])
histogram2(chromosome_vol_Symbiodinium_pilosum(:,1),chromosome_vol_Symbiodinium_pilosum(:,2),'FaceColor',[0.4940 0.1840 0.5560])
histogram2(chromosome_vol_Symbiodinium_minutum(:,1),chromosome_vol_Symbiodinium_minutum(:,2),'FaceColor',[0.9290 0.6940 0.1250])
histogram2(chromosome_vol_Brandtodinium_nutricula(:,1),chromosome_vol_Brandtodinium_nutricula(:,2),'FaceColor',[0.6350 0.0780 0.1840])
hold off
xlabel('Chromosome Volume [nm^3]')
set(gca,'YTickLabel',[]);
zlabel('Count')
ylabel('Cell')
ax = gca;
ax.FontSize = 24;
lgd=legend({'Symbiodinium microadriaticum','Symbiodinium pilosum','Symbiodinium minutum','Brandtodinium nutricula'});
legend boxoff
lgd.FontSize = 14;