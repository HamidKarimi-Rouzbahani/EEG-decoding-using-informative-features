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
'Cros Cor','Wavelet','Hilb Amp','Hilb Phs','Samples'};

series=2;  % 1:4
if series==1
    array=1:5;
    miny=-1.8;
elseif series==2
    array=6:14;
    miny=-2.68;
elseif series==3
    array=15:21;
    miny=-2.25;
elseif series==4
    array=22:26;
    miny=-1.8;
end
    
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 34];    

minx=-175;
maxx=975;
maxy=1;

accuracies_cued=nan*ones(26,231,4,10);
cued=1;
for Subject=[1:10]
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Cued_',num2str(cued),'_Band_',Bands{band},'_Wind_slid_Subject_',num2str(Subject),'.mat'],'accuracy');
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_cmplt.mat'],'accuracy');
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
    for feat=[1:26]
        Decoding_RT_correlaation_all(time,feat)=corr(nanmean(squeeze(nanmean(accuracies_cued(feat,time,:,:),3)),2),nanmean(Behavioural_RT_cued)','Type','Spearman');
    end
end
% save('Decoding_RT_correlaation_all_DS2.mat','Decoding_RT_correlaation_all')
% ccc
%% Plotting
figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(times,smooth(Decoding_RT_correlaation_all(:,Feature),20),'Color',colors{p},'linewidth',3);
    hold on;
end

line([minx maxx],[0 0],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');
%% Statistical analysis
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

% Decoding_Acc_correlaation_all_random=nan*ones(231,19,1000);
% Decoding_RT_correlaation_all_random=nan*ones(231,19,1000);
% for iteration=1:1000
%     for time=1:231
%         for feat=[1:26]
%             %         Decoding_Acc_correlaation_all_random(time,feat)=corr(nanmean([squeeze(nanmean(accuracies_cued(feat,time,:,:),3)) squeeze(nanmean(accuracies_uncued(feat,time,:,:),3))],2),nanmean([Behavioural_accuracy_cued;Behavioural_accuracy_uncued])','Type','Spearman');
%             Decoding_RT_correlaation_all_random(time,feat,iteration)=corr(nanmean(squeeze(nanmean(accuracies_cued(feat,time,:,:),3)),2),randsample(nanmean(Behavioural_RT_cued),10)','Type','Spearman');
%         end
%     end
% end
% save('random_RT_Behaav_correlations_cmplt.mat','Decoding_RT_correlaation_all_random')
load('random_RT_Behaav_correlations_cmplt.mat','Decoding_RT_correlaation_all_random');
threshold=0.9;
for feat=[1:26]
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


Baseline=-0.8;
steps=0.02;
distans=5; % times step
f=0;
Bayes=Bayes';
for feature=array
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==2
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{f},'linewidth',2,'markersize',5);
        elseif Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==4.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{f},'Color',colors{f},'linewidth',2,'markersize',5);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps+0.04;
    line([minx maxx],[baseline_temp baseline_temp]+1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+3.5*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-3.5*steps,'Color','k','linewidth',1);
end

if series==1
    
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
        {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)}},'EdgeColor','w','FontSize',14);

elseif series==2
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6),plot_line(7),plot_line(8),plot_line(9)],...
        {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)},Feat_names{array(6)},Feat_names{array(7)},Feat_names{array(8)},Feat_names{array(9)}},'EdgeColor','w','FontSize',10);

elseif series==3
     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6),plot_line(7)],...
        {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)},Feat_names{array(6)},Feat_names{array(7)}},'EdgeColor','w','FontSize',12);   
    
elseif series==4
    
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
        {Feat_names{array(1)},Feat_names{array(2)},Feat_names{array(3)},Feat_names{array(4)},Feat_names{array(5)}},'EdgeColor','w','FontSize',14);
end
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

data=nanmean(Decoding_RT_correlaation_all(35:end,:));

save('Mean_decoding_behav_correlation_each_feature_cmplt.mat','data');
ccc

%% Feature combinations

clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=2; % only DS2

series=1;  % 1:3
if series==1
    array=1:6;
    miny=-1.9;
elseif series==2
    array=[7:12];
    miny=-1.9;
elseif series==3
    array=13:17;
    miny=-1.7;
end
    
features=[1:9 12:19];    

minx=-175;
maxx=975;
maxy=1;

accuracies_cued=nan*ones(4,231,19,10);
Feat_Select_all=nan*ones(231,26,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
cued=1;
f=0;
for feat=features
    f=f+1;
    selection_method=listFS{f};
    for Subject=[1:10]
       load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy','Sel_feat');
        accuracy_orig=accuracy;
        
        try
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms.mat'],'accuracy');
            accuracy_tmp=smooth(nanmean(accuracy(1,:,206:end),2),4);
        catch
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy');
            accuracy_tmp=smooth(repmat(squeeze(nanmean(accuracy(1,:,47:53),2)),[4 1]),4);
            accuracy_tmp=accuracy_tmp(1:26);
        end
        
        accuracy=accuracy_orig;       
        % interpolation
        for cl=1:size(accuracy,2)
            accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
            accuracy_IP(1,cl,206:end)=accuracy_tmp;
            accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
        end
        
        accuracy=accuracy_IP;
        accuracies_cued(1,:,feat,Subject)=squeeze(nanmean(accuracy(1,[1 2 3],:),2));
        accuracies_cued(2,:,feat,Subject)=squeeze(nanmean(accuracy(1,[1 4 5],:),2));
        accuracies_cued(3,:,feat,Subject)=squeeze(nanmean(accuracy(1,[2 4 6],:),2));
        accuracies_cued(4,:,feat,Subject)=squeeze(nanmean(accuracy(1,[3 5 6],:),2));        
        [~,Behavioural_RT_cued(:,Subject)]= DatasetLoading_DS2_behaviour(Subject,cued); % category, percentage
    end
end

times=[-200:5:950]+25;
% Correlation across subjects
Decoding_RT_correlaation_all=nan*ones(231,19);
for time=1:length(times)
    f=0;
    for feat=features
        f=f+1;
        Decoding_RT_correlaation_all(time,f)=corr(nanmean(squeeze(nanmean(accuracies_cued(:,time,feat,:),1)),2),nanmean(Behavioural_RT_cued)','Type','Spearman');
    end
end
% save('Decoding_RT_correlaation_all__combinations_DS2.mat','Decoding_RT_correlaation_all')
figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(times,smooth(Decoding_RT_correlaation_all(:,Feature),20),'Color',colors{p},'linewidth',3);
    hold on;
end

line([minx maxx],[0 0],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');
% ccc
% % Statistical analysis
% Decoding_RT_correlaation_all_random=nan*ones(231,17,1000);
% for iteration=1:1000
%     for time=1:231
%         p=0;
%         for feat=features
%             p=p+1;
%             Decoding_RT_correlaation_all_random(time,p,iteration)=corr(nanmean(squeeze(nanmean(accuracies_cued(:,time,feat,:),1)),2),randsample(nanmean(Behavioural_RT_cued),10)','Type','Spearman');
%         end
%     end
%     iteration
% end
% save('random_RT_Behaav_correlations_feat_combin_cmplt.mat','Decoding_RT_correlaation_all_random')
% ccc

load('random_RT_Behaav_correlations_feat_combin_cmplt.mat','Decoding_RT_correlaation_all_random');
threshold=0.9;
for feat=1:17
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
    for e=1:17
        if Effects(Time,e)>0
            Bayes(Time,e)=4.5;
        elseif Effects(Time,e)<0
            Bayes(Time,e)=-0.5;
        elseif Effects(Time,e)==0
            Bayes(Time,e)=2;
        end
    end
end


Baseline=-0.8;
steps=0.02;
distans=4; % times step
f=0;
Bayes=Bayes';
for feature=array
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==2
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{f},'linewidth',2,'markersize',5);
        elseif Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==4.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{f},'Color',colors{f},'linewidth',2,'markersize',5);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps+0.04;
    line([minx maxx],[baseline_temp baseline_temp]+1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-1*steps,'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+3.5*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-3.5*steps,'Color','k','linewidth',1);
end

if series==1 || series==2
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)},listFS{array(6)}},'EdgeColor','w','FontSize',14);

elseif series==3
     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)}},'EdgeColor','w','FontSize',14);   
end
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



data=nanmean(Decoding_RT_correlaation_all(35:end,:));
save('Mean_decoding_behav_feat_combin_correlation_each_FS.mat','data');

