clc;
clear all;
% close all;

load('Max_decoding_accuracy_cmplt.mat');
data_Decod(:,1)=squeeze(nanmean(data(:,2,:),3))*100;
load('Mean_decoding_accuracy_cmplt.mat');
data_Decod(:,2)=squeeze(nanmean(data(:,2,:),3)*100);
load('First_above_chance_decoding_accuracy_cmplt.mat')
data_Decod(:,3)=squeeze(nanmean(data(:,2,:),3));
load('Max_time_decoding_accuracy_cmplt.mat');
data_Decod(:,4)=squeeze(nanmean(data(:,2,:),3));

load('Mean_decoding_behav_correlation_each_feature_cmplt.mat');
data_Corr=data';

% load('Mean_decoding_behav_feat_combned_correlation_each_feature_cmplt.mat');

Xdata={'Maximum Decoding Accuracy (%)','Average Decoding Accuracy (%)','Time of First Above-Chance Decoding (ms)','Time of Maximum Decoding (ms)'};
for i=1:4
    subplot_tmp=subplot(2,2,i);
    x=data_Decod(:,i);
    y=data_Corr;
    f=polyfit(x,y,1);
    scatter(x,y,80,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize', 11);

    y_est = polyval(f,x);
    hold on
    plot(x,y_est,'k--','LineWidth',3)
    hold off
    xlabel(Xdata{i})    
    ylabel('Average Correlation to Behavior (\rho)')';        
    set(subplot_tmp,'FontSize',14,'LineWidth',2);
    
    [r,p]=corr(data_Decod(:,i),data_Corr)
    text(max(data_Decod(:,i)),-0.2,['r = ',num2str(r,'%.2f')],'FontSize',16)
    text(max(data_Decod(:,i)),-0.25,['p = ',num2str(p,'%.2f')],'FontSize',16)
end
