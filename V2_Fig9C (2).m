clc;
clear all;
% close all;
load('Decoding_RT_correlaation_all_DS2.mat','Decoding_RT_correlaation_all')

% [a,numb]=min(min(Decoding_RT_correlaation_all(36:end,:)))
Decoding_RT_correlaation_all_both(1,:)=Decoding_RT_correlaation_all(:,23);
load('Decoding_RT_correlaation_all__combinations_DS2.mat','Decoding_RT_correlaation_all')
Decoding_RT_correlaation_all_both(2,:)=Decoding_RT_correlaation_all(:,8);
% [a,numb]=min(min(Decoding_RT_correlaation_all(36:end,:)))
% ccc
%% Plotting

figure;
colors={[0 0 0],[0 0.8 0]};
times=[-200:5:950]+25;

p=0;
for Feature=1:2
    p=p+1;
    plot_line{p}=plot(times,smooth(Decoding_RT_correlaation_all_both(p,:),20),'color',colors{p},'LineWidth',5);
    hold on;
end

minx=-175;
maxx=975;
maxy=1;
miny=-1;
line([minx maxx],[0 0],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');

legend([plot_line{1, 1},plot_line{1, 2}],{'Best individual feature (Wavelet)','Best combined FS method (laplacian)'},'EdgeColor','w','FontSize',20);
ylim([miny maxy])


ylabel(['Spearman''s Correlation to Behavior (\rho)'])
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',24,'LineWidth',4,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [-0.5 0 0.5 1],'YTickLabel',{'-0.5','0','0.5','1'},'XMinorTick','on');
xlim([minx maxx])
xtickangle(45);
