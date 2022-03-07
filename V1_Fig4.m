clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=1;

array=1:21;
features=[2:8 9 11:13 18 19 20 27 21:26];
chosen_features=features(array); %
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq'};

%% Significance
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
gca = axes('Position',[0.13 0.131 0.775 0.2]);
xtic=[1:5 7:15 17:23]*4;

for Dataset=1:3
    accuracies=nan*ones(35,231,10);
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,:,Subject)=nanmean(accuracy,2);
        for f=1:35
            accuracies(f,:,Subject)=smooth(accuracies(f,:,Subject),5);
        end
    end
  
    f=0;
    for feat=features
        f=f+1;
        significance(Dataset,f)=bf.ttest(nanmax(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
%           significance(Dataset,f)=bf.ttest(nanmean(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
    end
end

% Bayes stats againts chance
for Dataset=1:length(features)
    Effects=significance';
    for e=1:size(Effects,2)
        if Effects(Dataset,e)>10
            Bayes(Dataset,e)=2.5;
        elseif Effects(Dataset,e)>3 && Effects(Dataset,e)<=10
            Bayes(Dataset,e)=1.5;
        elseif Effects(Dataset,e)>1 && Effects(Dataset,e)<=3
            Bayes(Dataset,e)=0.5;
        elseif Effects(Dataset,e)<1 && Effects(Dataset,e)>=1/3
            Bayes(Dataset,e)=-0.5;
        elseif Effects(Dataset,e)<1/3 && Effects(Dataset,e)>=1/10
            Bayes(Dataset,e)=-1.5;
        elseif Effects(Dataset,e)<1/10
            Bayes(Dataset,e)=-2.5;
        end
    end
end
for Dataset=1:length(xtic)
    line([xtic(Dataset) xtic(Dataset)]+2,[0.39 0.43],'Color','k','linestyle',':','linewidth',1);
    hold on;
end

Baseline=0.45;
steps=0.005;
distans=1.5; % times step
for Dataset=1:size(Bayes,2)
    for f=1:size(Bayes,1)
        hold on;
        if Bayes(f,Dataset)==-0.5 || Bayes(f,Dataset)==0.5
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{9-Dataset},'linewidth',2,'markersize',7);
        elseif Bayes(Dataset,Dataset)~=0
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{9-Dataset},'Color',colors{9-Dataset},'linewidth',2,'markersize',7);
        end
    end
    baseline_temp=Baseline-(3*2+distans)*steps;
    line([-5 max(xtic)+5],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',1);
end

set(gca,'FontSize',18,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [1],'YTickLabel',{''});
ylabel({'Bayes';'Factor'})
xtickangle(45);
xlim([3 97])
ylim([0.39 0.43])
box off;


%% Bar plots
gca = axes('Position',[0.13 0.05 0.775 0.794]);
for Dataset=1:3
    accuracies=nan*ones(35,231,10);
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,:,Subject)=nanmean(accuracy,2);
        for f=1:35
            accuracies(f,:,Subject)=smooth(accuracies(f,:,Subject),5);
        end
    end
    f=0;
    for feat=features
        f=f+1;
        [data(f,Dataset,:),data_max(f,Dataset,:)]=nanmax(squeeze(accuracies(feat,31:end,:)));
%         data(f,Dataset,:)=nanmean(squeeze(accuracies(feat,31:end,:)));
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(data(f,Dataset,:)),'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(data(f,Dataset,:)),nanstd(data(f,Dataset,:))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
        significance(Dataset,f)=bf.ttest(squeeze(data(f,Dataset,:)),nanmean(squeeze(accuracies(feat,1:30,:)))');
    end
end

line([-5 max(xtic)+5],[0.5 0.5],'color','k','linestyle','--')
set(gca,'FontSize',18,'FontName','Calibri','XTick',xtic+2,'XTickLabel',{''},...
    'YTick',[0.5 0.55 0.6 0.65],'YTickLabel',{'0.5','0.55','0.6','0.65'});
xtickangle(45)
box off;
ylabel('Maximum Decoding Accuracy')
ylim([0.45 0.68]);
xlim([3 97])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);
% save('Mean_decoding_accuracy.mat','data')
% save('Max_decoding_accuracy.mat','data')
ccc

%% Significant times
clc;
clear all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=1;
features=[2:8 9 11:13 18 19 20 27 21:26];
xtic=[1:5 7:15 17:23]*4;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq'};

% ind_sig_total=nan*ones(21,3,10);
% for Dataset=1:3
%     accuracies=nan*ones(35,231,10);
%     if Dataset<3
%         accuracies_each_subj=nan*ones(35,6,231,10);
%     else
%         accuracies_each_subj=nan*ones(35,15,231,10);
%     end
%     for Subject=[1:10]
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
%         accuracies(:,:,Subject)=nanmean(accuracy,2);
%         accuracies_each_subj(:,:,:,Subject)=accuracy;
%         for f=1:35
%             accuracies(f,:,Subject)=smooth(accuracies(f,:,Subject),5);
%             for c=1:size(accuracy,2)
%                 accuracies_each_subj(f,c,:,Subject)=smooth(accuracies_each_subj(f,c,:,Subject),5);
%             end
%         end
%     end
%     
%     ff=0;
%     for feat=features
%         ff=ff+1;
%         [~,data_max(ff,Dataset,:)]=nanmax(squeeze(accuracies(feat,42:110,:)));
%     end
%     for Subject=[1:10]
%         significance=nan*ones(size(accuracies,1),size(accuracies,2));
%         for feature=1:size(accuracies,1)
%             for time=1:size(accuracies,2)
%                 significance(feature,time)=bf.ttest(squeeze(accuracies_each_subj(feature,:,time,Subject))',nanmean(squeeze(accuracies_each_subj(feature,:,1:30,Subject)),2));
%             end
%         end
%         Bayes=nan*ones(size(significance,1),size(significance,2));
%         for feature=1:size(significance,1)
%             Effects=significance;
%             for time=1:size(significance,2)
%                 if Effects(feature,time)>10
%                     Bayes(feature,time)=2.5;
%                 elseif Effects(feature,time)>3 && Effects(feature,time)<=10
%                     Bayes(feature,time)=1.5;
%                 elseif Effects(feature,time)>1 && Effects(feature,time)<=3
%                     Bayes(feature,time)=0.5;
%                 elseif Effects(feature,time)<1 && Effects(feature,time)>=1/3
%                     Bayes(feature,time)=-0.5;
%                 elseif Effects(feature,time)<1/3 && Effects(feature,time)>=1/10
%                     Bayes(feature,time)=-1.5;
%                 elseif Effects(feature,time)<1/10
%                     Bayes(feature,time)=-2.5;
%                 end
%                 Bayes_total(feature,time,Subject)=Bayes(feature,time);
%             end
%         end
%     end
%     
%     for Subject=[1:10]
%         f=0;
%         for feat=features
%             f=f+1;
%             tmp=find(Bayes_total(feat,42:100,Subject)>0.5);            
%             if length(tmp)>1
%                 for count=1:length(tmp)-1
%                     if tmp(count+1)-tmp(count)<2
%                         ind_sig(f,Dataset)=tmp(count);
%                         break;
%                     else
%                         ind_sig(f,Dataset)=nan;
%                     end
%                 end
%             else
%                 ind_sig(f,Dataset)=nan;
%             end
%             ind_sig_total(f,Dataset,Subject)=ind_sig(f,Dataset);
%         end
%     end
%     f=0;
%     for feat=features
%         f=f+1;
%         Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
%         hold on;
%         errorbar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,(nanstd(squeeze(ind_sig_total(f,Dataset,:)*5+30)'))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
%     end
% end
% save('First_above_chance_time.mat','ind_sig_total','data_max');

load('First_above_chance_time.mat');

for Dataset=1:3
    
    f=0;
    for feat=features
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,(nanstd(squeeze(ind_sig_total(f,Dataset,:)*5+30)'))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',18,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Time of First Above-Chance Decoding (ms)'})
xtickangle(45);
xlim([3 97])
box off;
ylim([0 700])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);

% data=ind_sig_total*5+30;
% save('First_above_chance_decoding_accuracy.mat','data')


% Peak times
figure;
gca = axes('Position',[0.13 0.05 0.775 0.794]);
for Dataset=1:3
    f=0;
    for feat=features
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),nanstd(data_max(f,Dataset,:)*5)./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',18,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Peak Decoding Time (ms)'})
xtickangle(45);
xlim([3 97])
ylim([0 700])
box off;
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);

% data=data_max*5+30;
% save('Max_time_decoding_accuracy.mat','data')

ccc
%% Cross-condition Significance Matrix 
data=ind_sig_total;
% data=data_max;
Dataset=3;

for feat1=1:21
    for feat2=1:21
        Significanc_mat(feat1,feat2)=bf.ttest(squeeze(data(feat1,Dataset,:)),squeeze(data(feat2,Dataset,:)));
        if feat1==feat2
            Significanc_mat(feat1,feat2)=1;
        end
    end
end
for feat1=1:21
    for feat2=1:21
        if Significanc_mat(feat1,feat2)>10
            Bayes_mat(feat1,feat2,Dataset)=6;
        elseif Significanc_mat(feat1,feat2)>3 && Significanc_mat(feat1,feat2)<=10
            Bayes_mat(feat1,feat2,Dataset)=5;
        elseif Significanc_mat(feat1,feat2)>1 && Significanc_mat(feat1,feat2)<=3
            Bayes_mat(feat1,feat2,Dataset)=4;
        elseif Significanc_mat(feat1,feat2)<1 && Significanc_mat(feat1,feat2)>=1/3
            Bayes_mat(feat1,feat2,Dataset)=3;
        elseif Significanc_mat(feat1,feat2)<1/3 && Significanc_mat(feat1,feat2)>=1/10
            Bayes_mat(feat1,feat2,Dataset)=2;
        elseif Significanc_mat(feat1,feat2)<1/10
            Bayes_mat(feat1,feat2,Dataset)=1;
        end
    end
end
subplot_tmp=subplot(1,1,1);
hold(subplot_tmp,'on');
image(Bayes_mat(:,:,Dataset),'Parent',subplot_tmp,'CDataMapping','scaled');
axis(subplot_tmp,'tight');
axis(subplot_tmp,'ij');
set(subplot_tmp,'CLim',[1 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');

xticks(1:21)
yticks(1:21)
xticklabels(Feat_names)
yticklabels(Feat_names)
ytickangle(45)
xtickangle(45)
colormap(parula(6));
c = colorbar ('YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
c.Ticks = [1+0.4 2.2 3.1 3.9 4.7 5.6]; %Create ticks
c.Label.String = 'Bayes Factors';
