clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=1;
MeanMax=1; %1=mean; 2=Max

array=1:28;
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 36 37 34];
chosen_features=features(array); %
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
    'Cross Corr','Wavelet','Hilb Amp','Hilb Phs','Amp Lock','Phs Lock','Orig Mag'};

%% Significance
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
gca = axes('Position',[0.16 0.35 0.85 0.1]);
xtic=[1:5 7:15 17:23 25:31]*4;

for Dataset=1:3
    accuracies=nan*ones(37,231,10);
    for Subject=[1:10]
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,:,Subject)=nanmean(accuracy,2);
        for f=1:37
            accuracies(f,:,Subject)=smooth(accuracies(f,:,Subject),5);
        end
    end
    
    f=0;
    for feat=features
        f=f+1;
        %% Do ou want mean of maximum of data
        if MeanMax==1
            significance(Dataset,f)=bf.ttest(nanmean(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
        else
            significance(Dataset,f)=bf.ttest(nanmax(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
        end
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
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{9-Dataset},'linewidth',2,'markersize',3);
        elseif Bayes(Dataset,Dataset)~=0
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{9-Dataset},'Color',colors{9-Dataset},'linewidth',2,'markersize',3);
        end
    end
    baseline_temp=Baseline-(3*2+distans)*steps;
    line([-5 max(xtic)+5],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',0.5);
end

set(gca,'FontSize',20,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [1],'YTickLabel',{''});
ylabel({'Bayes';'Factors'},'FontAngle','italic')
xtickangle(45);
xlim([3 130])
ylim([0.39 0.43])
box off;


%% Bar plots
gca = axes('Position',[0.16 0.45 0.85 0.55]);%x position, y position, width and height

for Dataset=1:size(Bayes,2)
    accuracies=nan*ones(37,231,10);
    for Subject=[1:10]
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,:,Subject)=nanmean(accuracy,2);
        for f=1:37
            accuracies(f,:,Subject)=smooth(accuracies(f,:,Subject),5);
        end
    end
    f=0;
    for feat=features
        f=f+1;
        %% do you want the max or mean parameters
        if MeanMax==1
            data(f,Dataset,:)=nanmean(squeeze(accuracies(feat,31:end,:)));
        else
            [data(f,Dataset,:),data_max(f,Dataset,:)]=nanmax(squeeze(accuracies(feat,31:end,:)));
        end
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(data(f,Dataset,:)),'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(data(f,Dataset,:)),nanstd(data(f,Dataset,:))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
        significance(Dataset,f)=bf.ttest(squeeze(data(f,Dataset,:)),nanmean(squeeze(accuracies(feat,1:30,:)))');
    end
end

line([-5 max(xtic)+5],[0.5 0.5],'color','k','linestyle','--')
set(gca,'FontSize',24,'FontName','Calibri','XTick',xtic+2,'XTickLabel',{''},...
    'YTick',[0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'});
xtickangle(45)
box off;
ylim([0.45 0.68]);
xlim([3 130])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','north','EdgeColor',[1 1 1]);
if MeanMax==1
    ylabel({'Average';'Decoding Accuracy (%)'})
    save('Revise_Mean_decoding_accuracy.mat','data')
else
    ylabel({'Maximum';'Decoding Accuracy (%)'})
    save('Revise_Max_decoding_accuracy.mat','data')
end

%% Significant times
clc;
clear all;
figure;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=2;
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 36 37 34];
xtic=[1:5 7:15 17:23 25:31]*4;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
    'Cross Corr','Wavelet','Hilb Amp','Hilb Phs','Amp Lock','Phs Lock','Orig Mag'};
gca = axes('Position',[0.16 0.45 0.85 0.55]);%x position, y position, width and height
 
% ind_sig_total=nan*ones(28,3,10);
% for Dataset=1:3
%     accuracies=nan*ones(37,231,10);
%     if Dataset<3
%         accuracies_each_subj=nan*ones(37,6,231,10);
%     else
%         accuracies_each_subj=nan*ones(37,15,231,10);
%     end
%     for Subject=[1:10]
%         load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
%         accuracies(:,:,Subject)=nanmean(accuracy,2);
%         accuracies_each_subj(:,:,:,Subject)=accuracy;
%         for f=1:37
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
% save('Revise_First_above_chance_time.mat','ind_sig_total','data_max');
% ccc
load('Revise_First_above_chance_time.mat');
for Dataset=1:3
    
    f=0;
    for feat=features
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,(nanstd(squeeze(ind_sig_total(f,Dataset,:)*5+30)'))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',20,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Time of First';'Above-Chance Decoding (ms)'})
xtickangle(45);
xlim([3 130])
box off;
ylim([0 350])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','north','EdgeColor',[1 1 1]);

data=ind_sig_total*5+30;
save('Revise_First_above_chance_decoding_accuracy.mat','data')


%% Peak times
figure;
gca = axes('Position',[0.16 0.45 0.85 0.55]);%x position, y position, width and height

for Dataset=1:3
    f=0;
    for feat=features
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),nanstd(data_max(f,Dataset,:)*5)./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',20,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',Feat_names,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Time of';'Maximum Decoding (ms)'})
xtickangle(45);
xlim([3 130])
ylim([0 350])
box off;
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','north','EdgeColor',[1 1 1]);


data=data_max*5+30;
save('Revise_Max_time_decoding_accuracy.mat','data')

ccc
%% Cross-condition Significance Matrix
clc;
clear all;
close all;
% load('Revise_Max_decoding_accuracy.mat','data')
% load('Revise_Mean_decoding_accuracy.mat','data')
% load('Revise_First_above_chance_decoding_accuracy.mat')
load('Revise_Max_time_decoding_accuracy.mat','data')

Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
    'Cross Corr','Wavelet','Hilb Amp','Hilb Phs','Amp Lock','Phs Lock','Orig Mag'};

for Dataset=1:3
    for feat1=1:28
        for feat2=1:28
            Significanc_mat(feat1,feat2)=bf.ttest(squeeze(data(feat1,Dataset,:)),squeeze(data(feat2,Dataset,:)));
            %         if feat1==feat2
            %             Significanc_mat(feat1,feat2)=1;
            %         end
            if feat1==feat2
                Significanc_mat(feat1,feat2)=0.01;
            end
            if feat1<feat2
                Significanc_mat(feat1,feat2)=0;
            end
            
        end
    end
    for feat1=1:28
        for feat2=1:28
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
            elseif Significanc_mat(feat1,feat2)<1/10 && Significanc_mat(feat1,feat2)~=0
                Bayes_mat(feat1,feat2,Dataset)=1;
            elseif Significanc_mat(feat1,feat2)==0
                Bayes_mat(feat1,feat2,Dataset)=0;
            end
        end
    end
    subplot_tmp=subplot(1,3,Dataset);
    hold(subplot_tmp,'on');
    image(Bayes_mat(:,:,Dataset),'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',[0 6],'DataAspectRatio',[1 1 1],'FontSize',7,'FontName','Calibri');
    
    xticks(1:28)
    yticks(1:28)
    xticklabels(Feat_names)
    yticklabels(Feat_names)
    ytickangle(45)
    xtickangle(45)
    colormap_mine=parula(6);
    colormap_mine=vertcat([1 1 1],colormap_mine);
    colormap(colormap_mine);
    title(['Dataset ',num2str(Dataset)],'FontSize',16)
end
% c = colorbar ('FontSize',20,'Limits',[0.9 6],'YTick',[1:6],'YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
% c.Ticks = [1+0.4 2.1 3 3.85 4.7 5.55]; %Create ticks
% c.Label.String = 'Bayes Factors';
%% Correlation between parameters
clc;
clear all;
% close all;
for Dataset=1:3
    load('Revise_Max_decoding_accuracy.mat');
    data_Decod(:,1)=squeeze(nanmean(data(:,Dataset,:),3))*100;
    load('Revise_Mean_decoding_accuracy.mat');
    data_Decod(:,2)=squeeze(nanmean(data(:,Dataset,:),3)*100);
    load('Revise_First_above_chance_decoding_accuracy.mat')
    data_Decod(:,3)=squeeze(nanmean(data(:,Dataset,:),3));
    load('Revise_Max_time_decoding_accuracy.mat');
    data_Decod(:,4)=squeeze(nanmean(data(:,Dataset,:),3));
    
    
    Xdata={'Maximum Decoding Accuracy (%)','Average Decoding Accuracy (%)'};
    subplot_tmp=subplot(2,3,Dataset);
    x=data_Decod(:,2);
    y=data_Decod(:,1);
    f=polyfit(x,y,1);
    scatter(x,y,80,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize', 11);
    
    y_est = polyval(f,x);
    hold on
    plot(x,y_est,'k--','LineWidth',3)
    hold off
    xlabel(Xdata{2})
    ylabel(Xdata{1});
    set(subplot_tmp,'FontSize',14,'LineWidth',2);
    
    [r,p]=corr(data_Decod(:,2),data_Decod(:,1))
    text(52.5,55,['r = ',num2str(r,'%.2f')],'FontSize',20)
    text(52.5,52.5,['p = ',num2str(p,'%.2f')],'FontSize',20)
    xlim([48 56])
    ylim([50 70])
    title(['Dataset ',num2str(Dataset)])
    
    Xdata={'Time of First Above-Chance Decoding (ms)','Time of Maximum Decoding (ms)'};
    subplot_tmp=subplot(2,3,Dataset+3);
    x=data_Decod(:,3);
    y=data_Decod(:,4);
    f=polyfit(x,y,1);
    scatter(x,y,80,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize', 11);
    
    y_est = polyval(f,x);
    hold on
    plot(x,y_est,'k--','LineWidth',3)
    hold off
    xlabel(Xdata{1})
    ylabel(Xdata{2});
    set(subplot_tmp,'FontSize',12,'LineWidth',2);
    
    [r,p]=corr(data_Decod(:,3),data_Decod(:,4));
    text(150,130,['r = ',num2str(r,'%.2f')],'FontSize',20)
    text(150,115,['p = ',num2str(p,'%.2f')],'FontSize',20)
    xlim([50 200])
    ylim([100 250])
end