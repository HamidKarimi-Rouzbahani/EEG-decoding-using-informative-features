%% 1=Features one by one
clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=2;
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
'Cross Corr','Wavelet','Hilb Amp','Hilb Phs','Amp Lock','Phs Lock','Orig Mag'};

series=4;  % 1:4
if series==1
    array=1:5;
    miny=-1.97;
elseif series==2
    array=6:14;
    miny=-2.85;
elseif series==3
    array=15:21;
    miny=-2.4;
elseif series==4
    array=22:28;
    miny=-2.4;
end
    
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 36 37 34];    

minx=-175;
maxx=975;
maxy=1;

accuracies_cued=nan*ones(28,231,4,10);
cued=1;
for Subject=[1:10]
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Cued_',num2str(cued),'_Band_',Bands{band},'_Wind_slid_Subject_',num2str(Subject),'.mat'],'accuracy');
    load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracy(accuracy==0)=nan;
    accuracies_cued(:,:,1,Subject)=squeeze(nanmean(accuracy(features,[1 2 3],:),2));
    accuracies_cued(:,:,2,Subject)=squeeze(nanmean(accuracy(features,[1 4 5],:),2));
    accuracies_cued(:,:,3,Subject)=squeeze(nanmean(accuracy(features,[2 4 6],:),2));
    accuracies_cued(:,:,4,Subject)=squeeze(nanmean(accuracy(features,[3 5 6],:),2));
    [~,Behavioural_RT_cued(:,Subject)]= DatasetLoading_DS2_behaviour(Subject,cued); % category, percentage
end
% Correlation across subjects
Decoding_RT_correlaation_all=nan*ones(231,26);
for time=1:231
    for feat=[1:28]
        Decoding_RT_correlaation_all(time,feat)=corr(nanmean(squeeze(nanmean(accuracies_cued(feat,time,:,:),3)),2),nanmean(Behavioural_RT_cued)','Type','Spearman');
    end
end
% save('Revise_decoding_RT_correlaation_all_DS2.mat','Decoding_RT_correlaation_all')
% ccc
%% Plotting
figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(times,smooth(Decoding_RT_correlaation_all(:,Feature),20),'Color',colors{p},'linewidth',2);
    hold on;
end

line([minx maxx],[0 0],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');
%% Statistical analysis
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

% Decoding_Acc_correlaation_all_random=nan*ones(231,28,1000);
% Decoding_RT_correlaation_all_random=nan*ones(231,28,1000);
% for iteration=1:1000
%     for time=1:231
%         for feat=[1:28]
%             %         Decoding_Acc_correlaation_all_random(time,feat)=corr(nanmean([squeeze(nanmean(accuracies_cued(feat,time,:,:),3)) squeeze(nanmean(accuracies_uncued(feat,time,:,:),3))],2),nanmean([Behavioural_accuracy_cued;Behavioural_accuracy_uncued])','Type','Spearman');
%             Decoding_RT_correlaation_all_random(time,feat,iteration)=corr(nanmean(squeeze(nanmean(accuracies_cued(feat,time,:,:),3)),2),randsample(nanmean(Behavioural_RT_cued),10)','Type','Spearman');
%         end
%     end
%     iteration
% end
% save('Revise_random_RT_Behaav_correlations_cmplt.mat','Decoding_RT_correlaation_all_random')
load('Revise_random_RT_Behaav_correlations_cmplt.mat','Decoding_RT_correlaation_all_random');
threshold=0.9;
for feat=[1:28]
    for time=1:231      
        if Decoding_RT_correlaation_all(time,feat)>=0 && sum(Decoding_RT_correlaation_all(time,feat)>Decoding_RT_correlaation_all_random(time,feat,:))>threshold*size(Decoding_RT_correlaation_all_random,3)
            significance(time,feat)=20;
        elseif Decoding_RT_correlaation_all(time,feat)<0 && sum(Decoding_RT_correlaation_all(time,feat)<Decoding_RT_correlaation_all_random(time,feat,:))>threshold*size(Decoding_RT_correlaation_all_random,3)
            significance(time,feat)=-20;
        else
            significance(time,feat)=0;
        end
    end
end
for Time=1:length(significance)
    Effects=significance;
    for e=1:size(Effects,2)
        if Effects(Time,e)>0
            Bayes(Time,e)=4.5;
        elseif Effects(Time,e)<0
            Bayes(Time,e)=-0.5;
        elseif Effects(Time,e)==0
            Bayes(Time,e)=2;
        end
    end
end


Baseline=-0.9;
steps=0.02;
distans=5; % times step
f=0;
Bayes=Bayes';
for feature=array
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==2
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{f},'linewidth',2,'markersize',4);
        elseif Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==4.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{f},'Color',colors{f},'linewidth',2,'markersize',4);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps+0.04;
    line([minx maxx],[baseline_temp baseline_temp]+1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+3.5*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-3.5*steps,'Color','k','linewidth',1);
end

% if series==1
%     
%     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
%         {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)}},'EdgeColor','w','FontSize',14);
% 
% elseif series==2
%     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6),plot_line(7),plot_line(8),plot_line(9)],...
%         {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)},Feat_names{array(6)},Feat_names{array(7)},Feat_names{array(8)},Feat_names{array(9)}},'EdgeColor','w','FontSize',10);
% 
% elseif series==3
%      legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6),plot_line(7)],...
%         {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)},Feat_names{array(6)},Feat_names{array(7)}},'EdgeColor','w','FontSize',12);   
%     
% elseif series==4
%     
%      legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6),plot_line(7)],...
%         {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)},Feat_names{array(6)},Feat_names{array(7)}},'EdgeColor','w','FontSize',12);   
%     
% end
ylim([miny maxy])

ylabel(['Spearman''s Correlation to Behavior (\rho)'])
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',19,'LineWidth',3,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [-0.5 0 0.5 1],'YTickLabel',{'-0.5','0','0.5','1'},'XMinorTick','on');
xlim([minx maxx])
xtickangle(45);

data=nanmean(Decoding_RT_correlaation_all(35:end,:));

save('Revise_mean_decoding_behav_correlation_each_feature_cmplt.mat','data');
ccc
%%
clc;
clear all;
% close all;

load('Revise_Max_decoding_accuracy.mat');
data_Decod(:,1)=squeeze(nanmean(data(:,2,:),3))*100;
load('Revise_Mean_decoding_accuracy.mat');
data_Decod(:,2)=squeeze(nanmean(data(:,2,:),3)*100);
load('Revise_First_above_chance_decoding_accuracy.mat')
data_Decod(:,3)=squeeze(nanmean(data(:,2,:),3));
load('Revise_Max_time_decoding_accuracy.mat');
data_Decod(:,4)=squeeze(nanmean(data(:,2,:),3));


load('Revise_mean_decoding_behav_correlation_each_feature_cmplt.mat');
data_Corr=data';

% load('Mean_decoding_behav_feat_combned_correlation_each_feature.mat');

Xdata={'Maximum Decoding Accuracy (%)','Average Decoding Accuracy (%)','Time of First Above-Chance Decoding (ms)','Time of Maximum Decoding (ms)'};
for i=1:4
    subplot_tmp=subplot(2,2,i);
    x=data_Decod(:,i);
    y=data_Corr;
    f=polyfit(x,y,1);
    scatter(x,y,50,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize', 11);
    
    y_est = polyval(f,x);
    hold on
    plot(x,y_est,'k--','LineWidth',3)
    hold off
    xlabel(Xdata{i})
    ylabel('Average Correlation to Behavior (\rho)')';
    set(subplot_tmp,'FontSize',11,'LineWidth',2);
    
    [r,p]=corr(data_Decod(:,i),data_Corr)
    text(max(data_Decod(:,i)),-0.2,['r = ',num2str(r,'%.2f')],'FontSize',16)
    text(max(data_Decod(:,i)),-0.25,['p = ',num2str(p,'%.2f')],'FontSize',16)
end
