clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=1;

minx=-150;
maxx=950;
maxy=0.665;
miny=0.395;


accuracies1=nan*ones(231,10);
times1=[-200:5:950]+2.5;
accuracies2=nan*ones(231,10);
times2=[-200:5:950]+25;
accuracies3=nan*ones(231,10);
times3=[-200:5:950]+50;
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_5msWind.mat'],'accuracy');
    accuracies1(:,Subject)=nanmean(accuracy(2,:,:),2);
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies2(:,Subject)=nanmean(accuracy(2,:,:),2);
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_100msWind.mat'],'accuracy');
    %     accuracies3(:,Subject)=[squeeze(nanmean(accuracy(2,:,:),2));nan(10,1)];
    accuracies3(:,Subject)=[squeeze(nanmean(accuracy(2,:,:),2));squeeze(nanmean(accuracy(2,:,212:221),2))];
end

colors={[0.8 0 0],[0 0 0],[0 0 0.8]};
figure;
times=[-175:5:975];
for subj=1:10
    accuraciesT(1,:,subj)=interp1(times1,accuracies1(:,subj),times,'spline');
    accuraciesT(2,:,subj)=interp1(times2,accuracies2(:,subj),times,'spline');
    accuraciesT(3,:,subj)=interp1(times3,accuracies3(:,subj),times,'spline');
end

plot_line{1}=shadedErrorBar(times,smooth(nanmean(accuraciesT(1,:,:),3),5),smooth(nanstd(squeeze(accuraciesT(1,:,:))'),5)./sqrt(10),{'color',colors{1},'LineWidth',3},1);
hold on;
plot_line{2}=shadedErrorBar(times,smooth(nanmean(accuraciesT(2,:,:),3),5),smooth(nanstd(squeeze(accuraciesT(2,:,:))'),5)./sqrt(10),{'color',colors{2},'LineWidth',3},1);
plot_line{3}=shadedErrorBar(times,smooth(nanmean(accuraciesT(3,:,:),3),5),smooth(nanstd(squeeze(accuraciesT(3,:,:))'),5)./sqrt(10),{'color',colors{3},'LineWidth',3},1);
line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');

%% sigfnificance
for time=1:231
    significance(1,time)=bf.ttest(squeeze(accuraciesT(1,time,:)),squeeze(accuraciesT(2,time,:)));
    significance(2,time)=bf.ttest(squeeze(accuraciesT(3,time,:)),squeeze(accuraciesT(2,time,:)));
end

% Bayes stats againts chance
for feature=1:size(significance,1)
    Effects=significance;
    for time=1:size(significance,2)
        if Effects(feature,time)>10
            Bayes(feature,time)=2.5;
        elseif Effects(feature,time)>3 && Effects(feature,time)<=10
            Bayes(feature,time)=1.5;
        elseif Effects(feature,time)>1 && Effects(feature,time)<=3
            Bayes(feature,time)=0.5;
        elseif Effects(feature,time)<1 && Effects(feature,time)>=1/3
            Bayes(feature,time)=-0.5;
        elseif Effects(feature,time)<1/3 && Effects(feature,time)>=1/10
            Bayes(feature,time)=-1.5;
        elseif Effects(feature,time)<1/10
            Bayes(feature,time)=-2.5;
        end
    end
end

times=[-200:5:950]+50;
Baseline=0.47;
steps=0.005;
distans=2; % times step
f=0;
colors={[0.8 0 0],[0 0 0.8]};
for feature=1:2
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==0.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{feature},'linewidth',2,'markersize',4);
        elseif Bayes(feature,windoww)~=0
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{feature},'Color',colors{feature},'linewidth',2,'markersize',4);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps;
    line([minx maxx],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',1);
end


legend([plot_line{1, 1}.mainLine,plot_line{1, 2}.mainLine,plot_line{1, 3}.mainLine],{'Window length = 5 ms','Window length = 50 ms','Window length = 100 ms'},'EdgeColor','w','FontSize',20);

ylim([miny maxy])

ylabel('Decoding Accuracy (%)')
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',24,'LineWidth',4,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'},'XMinorTick','on');
xtickangle(45);
xlim([minx maxx])


