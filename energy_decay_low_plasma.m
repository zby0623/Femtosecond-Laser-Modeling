close all;
NP=[5 10 20 50 100];
legend_string=string(zeros(1,length(NP)));
for i=1:length(NP)
    load("filament_start_pulse50fs_0326_lowplasma_P_Pcr_"+num2str(NP(i))+".mat");
    plot((0:1*dz:z)/filament_start_loc1(1),Energy(1:1:round(z/dz+1))/max(Energy),'LineWidth',2);
    hold on;
    legend_string(i)="P_{in}/P_{cr}="+num2str(NP(i));
end
xlabel('Z/L_{c}');
ylabel('Energy/Original Energy');
set(gca,'FontSize',20);
set(gca,'FontWeight','Bold')
legend(legend_string,'Location','Southwest');
axis([0 6 0.6 1]);