clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Datasets_used=1:3;
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
method=[1:9 12:19];    

%% Significance
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
gca = axes('Position',[0.13 0.131 0.775 0.2]);
xtic=[1:17]*4;

for Dataset=Datasets_used
    accuracies=nan*ones(19,231,10);
    f=0;
    for feat=method
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
        if strcmp(selection_method,'fisher')
            accuracy=accuracy-repmat(nanmean(nanmean(accuracy(1,:,1:10),2),3),[size(accuracy)])+0.5;
        end
        
        % interpolation
            for cl=1:size(accuracy,2)
                accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
                accuracy_IP(1,cl,206:end)=accuracy_tmp;
                accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
            end
            accuracy=accuracy_IP;

            accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));
            for g=1:length(method)
                accuracies(g,:,Subject)=smooth(accuracies(g,:,Subject),5);
            end
        end
    end
    
    f=0;
    for feat=1:length(method)
        f=f+1;
%         significance(Dataset,f)=bf.ttest(nanmax(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
        significance(Dataset,f)=bf.ttest(nanmean(squeeze(accuracies(feat,31:end,:)))',nanmean(squeeze(accuracies(feat,1:30,:)))');
    end
end

% Bayes stats againts chance
for Dataset=1:length(method)
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
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{9-Dataset},'linewidth',2,'markersize',4);
        elseif Bayes(Dataset,Dataset)~=0
            plots(Dataset)=plot(xtic(f)+Dataset,Bayes(f,Dataset).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{9-Dataset},'Color',colors{9-Dataset},'linewidth',2,'markersize',4);
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
    [xtic]+2,'XTickLabel',listFS,'YTick',...
    [1],'YTickLabel',{''});
ylabel({'Bayes';'Factors'})
xtickangle(45);
xlim([3 73])
ylim([0.39 0.43])
box off;


%% Bar plots
gca = axes('Position',[0.13 0.05 0.775 0.794]);
for Dataset=Datasets_used
    accuracies=nan*ones(19,231,10);
    f=0;
    for feat=method
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
        if strcmp(selection_method,'fisher')
            accuracy=accuracy-repmat(nanmean(nanmean(accuracy(1,:,1:10),2),3),[size(accuracy)])+0.5;
        end
            % interpolation
            for cl=1:size(accuracy,2)
                accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
            accuracy_IP(1,cl,206:end)=accuracy_tmp;
                accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
            end
            accuracy=accuracy_IP;
            
            accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));
            for g=1:length(method)
                accuracies(g,:,Subject)=smooth(accuracies(g,:,Subject),5);
            end
        end
    end
    f=0;
    for feat=1:length(method)
        f=f+1;
%         [data(f,Dataset,:),data_max(f,Dataset,:)]=nanmax(squeeze(accuracies(feat,31:end,:)));
        data(f,Dataset,:)=nanmean(squeeze(accuracies(feat,31:end,:)));

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
ylabel('Average Decoding Accuracy (%)')
% ylabel('Maximum Decoding Accuracy (%)')
ylim([0.45 0.68]);
xlim([3 73])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);
save('Mean_decoding_accuracy_feat_combined_cmplt.mat','data')
% save('Max_decoding_accuracy_feat_combined_cmplt.mat','data')

ccc

%% Significant times
clc;
clear all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
xtic=[1:17]*4;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
Datasets_used=1:3;
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
method=[1:9 12:19];    
% 
% ind_sig_total=nan*ones(17,3,10);
% for Dataset=Datasets_used
%     accuracies=nan*ones(17,231,10);
%     if Dataset<3
%         accuracies_each_subj=nan*ones(17,6,231,10);
%     else
%         accuracies_each_subj=nan*ones(17,15,231,10);
%     end
%     f=0;
%     for feat=method
%         f=f+1;
%         selection_method=listFS{f};
%         for Subject=[1:10]
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy','Sel_feat');
%         accuracy_orig=accuracy;
%         
%         try
%             load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms.mat'],'accuracy');
%             accuracy_tmp=smooth(nanmean(accuracy(1,:,206:end),2),4);
%         catch
%             load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy');
%             accuracy_tmp=smooth(repmat(squeeze(nanmean(accuracy(1,:,47:53),2)),[4 1]),4);
%             accuracy_tmp=accuracy_tmp(1:26);
%         end
%         
%         accuracy=accuracy_orig;
%         if strcmp(selection_method,'fisher')
%             accuracy=accuracy-repmat(nanmean(nanmean(accuracy(1,:,1:10),2),3),[size(accuracy)])+0.5;
%         end
%         
%         % interpolation
%             for cl=1:size(accuracy,2)
%                 accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
%                 accuracy_IP(1,cl,206:end)=accuracy_tmp;
%                 accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
%             end
%             accuracy=accuracy_IP;
%             accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));            
%             accuracies_each_subj(f,:,:,Subject)=squeeze(accuracy);
%             for g=1:17
%                 accuracies(g,:,Subject)=smooth(accuracies(g,:,Subject),5);
%                 for c=1:size(accuracy,2)
%                     accuracies_each_subj(g,c,:,Subject)=smooth(accuracies_each_subj(g,c,:,Subject),5);
%                 end
%             end
%         end
%     end
%     ff=0;
%     for feat=1:length(method)
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
%     for Subject=[1:10]
%         f=0;
%         for feat=1:length(method)
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
%     
%     f=0;
%     for feat=1:length(method)
%         f=f+1;
%         Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
%         hold on;
%         errorbar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,(nanstd(squeeze(ind_sig_total(f,Dataset,:)*5+30)'))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
%     end
% end
% ccc
% save('First_above_chance_time_feat_combined_cmplt.mat','ind_sig_total','data_max');

load('First_above_chance_time_feat_combined_cmplt.mat');
for Dataset=Datasets_used
    
    f=0;
    for feat=1:length(method)
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(ind_sig_total(f,Dataset,:),3)*5+30,(nanstd(squeeze(ind_sig_total(f,Dataset,:)*5+30)'))./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',24,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',listFS,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Time of First Above-Chance Decoding (ms)'})
xtickangle(45);
xlim([3 73])
box off;
ylim([0 350])
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);
% legend([Bars(1),Bars(2)],{'Dataset 1','Dataset 2'},'location','northeast','EdgeColor',[1 1 1]);

data=ind_sig_total*5+30;
save('First_above_chance_decoding_accuracy_feat_combined_cmplt.mat','data')
ccc

% Peak times
figure;
gca = axes('Position',[0.13 0.05 0.775 0.794]);
for Dataset=1:3
    f=0;
    for feat=1:length(method)
        f=f+1;
        Bars(Dataset)=bar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),'facecolor',colors{9-Dataset},'edgecolor','none','LineWidth',0.1);
        hold on;
        errorbar(xtic(f)+Dataset,nanmean(data_max(f,Dataset,:)*5+30),nanstd(data_max(f,Dataset,:)*5)./sqrt(10),'linewidth',2,'color',colors{9-Dataset},'CapSize',0,'LineStyle','none')
    end
end
set(gca,'FontSize',24,'FontName','Calibri','XTick',...
    [xtic]+2,'XTickLabel',listFS,'YTick',...
    [100:100:800],'YTickLabel',{'100','200','300','400','500','600','700','800'});
ylabel({'Time of Maximum Decoding (ms)'})
xtickangle(45);
xlim([3 73])
ylim([0 350])
box off;
legend([Bars(1),Bars(2),Bars(3)],{'Dataset 1','Dataset 2','Dataset 3'},'location','northeast','EdgeColor',[1 1 1]);
% legend([Bars(1),Bars(2)],{'Dataset 1','Dataset 2'},'location','northeast','EdgeColor',[1 1 1]);

data=data_max*5+30;
save('Max_time_decoding_accuracy_feat_combined_cmplt.mat','data')

%% Cross-condition Significance Matrix 
clc;
clear all;
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
% load('Max_decoding_accuracy_feat_combined_cmplt.mat','data')
% load('Mean_decoding_accuracy_feat_combined_cmplt.mat','data')
% load('First_above_chance_decoding_accuracy_feat_combined_cmplt.mat','data')
load('Max_time_decoding_accuracy_feat_combined_cmplt.mat','data')
for DS=1:3
    [A,FS]=sort(nanmean(data(:,DS,:),3));
    listFS{FS(1:3)} % min
%     listFS{FS(end-2:end)} % max

end
ccc

method=[1:9 12:19];    
figure;

for Dataset=1:3
    for feat1=1:length(method)
        for feat2=1:length(method)
            Significanc_mat(feat1,feat2)=bf.ttest(squeeze(data(feat1,Dataset,:)),squeeze(data(feat2,Dataset,:)));
            if feat1==feat2
                Significanc_mat(feat1,feat2)=0.01;
            end
            if feat1<feat2
                Significanc_mat(feat1,feat2)=0;
            end         
            if feat1>feat2 && isnan(Significanc_mat(feat1,feat2))
                Significanc_mat(feat1,feat2)=0.2;
            end
        end
    end
    for feat1=1:length(method)
        for feat2=1:length(method)
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
    set(subplot_tmp,'CLim',[0 6],'DataAspectRatio',[1 1 1],'FontSize',14,'FontName','Calibri');
    
    xticks(1:17)
    yticks(1:17)
    xticklabels(listFS)
    yticklabels(listFS)
    ytickangle(45)
    xtickangle(45)
    colormap_mine=parula(6);
    colormap_mine=vertcat([1 1 1],colormap_mine);
    colormap(colormap_mine);
    title(['Dataset ',num2str(Dataset)],'FontSize',16)
end

%% Correlation between parameters
clc;
clear all;
% close all;

for Dataset=1:3

load('Max_decoding_accuracy_feat_combined_cmplt.mat');
data_Decod(:,1)=squeeze(nanmean(data(:,Dataset,:),3))*100;
load('Mean_decoding_accuracy_feat_combined_cmplt.mat');
data_Decod(:,2)=squeeze(nanmean(data(:,Dataset,:),3)*100);
load('First_above_chance_decoding_accuracy_feat_combined_cmplt.mat')
data_Decod(:,3)=squeeze(nanmean(data(:,Dataset,:),3));
load('Max_time_decoding_accuracy_feat_combined_cmplt.mat');
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
text(52,53,['r = ',num2str(r,'%.2f')],'FontSize',20)
text(52,52.3,['p = ',num2str(p,'%.2f')],'FontSize',20)
xlim([50 53])
ylim([51 60])
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
set(subplot_tmp,'FontSize',14,'LineWidth',2);

[r,p]=corr(data_Decod(:,3),data_Decod(:,4))
text(160,230,['r = ',num2str(r,'%.2f')],'FontSize',20)
text(160,220,['p = ',num2str(p,'%.2f')],'FontSize',20)
xlim([70 200])
ylim([140 250])
end