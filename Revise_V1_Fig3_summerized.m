clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;

Dataset=3;
array=[1 23 28];
miny=0.24;

all_features_sorted=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 36 37 34];
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
    'Cross Corr','Wavelet','Hilb Amp','Hilb Phs','Amp Lock','Phs Lock','Orig Mag'};
chosen_features=all_features_sorted(array); %

minx=-175;
maxx=975;
maxy=0.66;

accuracies=nan*ones(37,231,10);
for Subject=[1:10]
    load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies(:,:,Subject)=nanmean(accuracy,2);
end
%% Plotting

figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

p=0;
for Feature=chosen_features
    p=p+1;
    plot_line{p}=shadedErrorBar(times,nanmean(squeeze((accuracies(Feature,:,:))),2),nanstd(squeeze((accuracies(Feature,:,:)))')'./sqrt(10),{'color',colors{p},'LineWidth',2},1);
    hold on;
end

line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');


%% sigfnificance
f=0;
for feature=chosen_features
    f=f+1;
    for time=1:size(accuracies,2)
        significance(f,time)=bf.ttest(squeeze(accuracies(feature,time,:)),squeeze(nanmean(accuracies(feature,1:30,:),2)));
    end
end

combins=nchoosek([1:length(array)],2);
for cobs=1:size(combins,1)
    f=f+1;
    for time=1:size(accuracies,2)
        significance(f,time)=bf.ttest(squeeze(accuracies(chosen_features(combins(cobs,1)),time,:)),squeeze(accuracies(chosen_features(combins(cobs,2)),time,:)));
    end
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


Baseline=0.47;
steps=0.005;
distans=2; % times step
f=0;
for feature=1:size(Bayes,1)
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(f,windoww)==-0.5 || Bayes(f,windoww)==0.5
            plots(f)=plot(times(windoww),Bayes(f,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{f},'linewidth',2,'markersize',2);
        elseif Bayes(f,windoww)~=0
            plots(f)=plot(times(windoww),Bayes(f,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{f},'Color',colors{f},'linewidth',2,'markersize',2);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps;
    line([minx maxx],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',0.5);
    line([minx maxx],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',0.5);
end

% legend([plot_line{1,1}.mainLine,plot_line{1,2}.mainLine,plot_line{1,3}.mainLine],...
%     {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)}},'EdgeColor','w','FontSize',14);

ylim([miny maxy])

ylabel('Decoding Accuracy (%)')
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',19,'LineWidth',3,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'},'XMinorTick','on');
xtickangle(45);
xlim([minx maxx])

