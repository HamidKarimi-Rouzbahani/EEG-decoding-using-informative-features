%% Feature combinations

clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=2; % only DS2

series=3;  % 1:3


if series==1
    array=1:6;
elseif series==2
    array=[7:12];
elseif series==3
    array=13:17;
end

features=[1:9 12:19];

minx=-175;
maxx=975;
maxy=.75;
miny=-.75;

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
        ccc
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
gca = axes('Position',[0.22 0.3 0.775 0.65]);
colors={[0 0.8 0.8],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
smoothing_rate=20;
p=0;
for Feature=array
    p=p+1;
    smoothened_data(p,:)=smooth(Decoding_RT_correlaation_all(:,Feature),smoothing_rate);
%     smoothened_data(p,:)=smoothdata(Decoding_RT_correlaation_all(:,Feature));
    plot(times,smoothened_data(p,:),'Color',colors{p},'linewidth',1);
    hold on;
end

line([minx maxx],[0 0],'LineWidth',1,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1,'Color','k','LineStyle','--');

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
% % ccc

load('random_RT_Behaav_correlations_feat_combin_cmplt.mat','Decoding_RT_correlaation_all_random');
threshold=0.9;
p=0;
for feat=array
    p=p+1;
    for time=1:231
        if Decoding_RT_correlaation_all(time,feat)>=0 && sum(Decoding_RT_correlaation_all(time,feat)>Decoding_RT_correlaation_all_random(time,feat,:))>threshold*size(Decoding_RT_correlaation_all_random,3)
            significance(time,p)=smoothened_data(p,time);
        elseif Decoding_RT_correlaation_all(time,feat)<0 && sum(Decoding_RT_correlaation_all(time,feat)<Decoding_RT_correlaation_all_random(time,feat,:))>threshold*size(Decoding_RT_correlaation_all_random,3)
            significance(time,p)=smoothened_data(p,time);
        else
            significance(time,p)=nan;
        end
    end
end

p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(times,significance(:,p),'Color',colors{p},'linewidth',3);
    hold on;
end

% if series==1 || series==2
%     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6)],...
%         {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)},listFS{array(6)}},'EdgeColor','none','FontSize',10,'location','north');
%     
% elseif series==3
%     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
%         {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)}},'EdgeColor','none','FontSize',10,'location','southeast');
% end
ylim([miny maxy])

ylabel({'Spearman''s Correlation';'to Behavior (\rho)'})
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',16,'LineWidth',1,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [-0.75 -0.375 0 0.375 0.75],'ycolor','k','xcolor','k','YTickLabel',{'-.75','-0.375','0','0.375','0.75'},'XMinorTick','on');
xlim([minx maxx])
xtickangle(45);
% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
pdf_paper_size         = [20 20];
fig.PaperSize       = pdf_paper_size;
pdf_print_resolution   = '-r300';
pdf_file_name=['Comb_paper_Fig4A_Series_',num2str(series)];
print(['Z:\projects\Hamid\Projects\Mojgan\Analyses\Combined_features\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)
