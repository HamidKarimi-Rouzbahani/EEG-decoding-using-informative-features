clc;
clear all;
close all;

load('Max_decoding_accuracy.mat');
data_Decod(:,1)=squeeze(nanmean(data(:,2,:),3));
load('Mean_decoding_accuracy.mat');
data_Decod(:,2)=squeeze(nanmean(data(:,2,:),3));
load('First_above_chance_decoding_accuracy.mat')
data_Decod(:,3)=squeeze(nanmean(data(:,2,:),3));
load('Max_time_decoding_accuracy.mat')
data_Decod(:,4)=squeeze(nanmean(data(:,2,:),3));

load('Mean_decoding_behav_correlation_each_feature.mat');
data_Corr=data';

% load('Mean_decoding_behav_feat_combned_correlation_each_feature.mat');

Xdata={'Maximum Decoding Accuracy','Average Decoding Accuracy','Time to First Decoding','Time to Maximum Decoding'};
for i=1:4
    subplot(2,2,i)
    x=data_Decod(:,i);
    y=data_Corr;
    f=polyfit(x,y,1);
    scatter(x,y,'marker','o','MarkerFaceColor','k');
    set(gca,'fontsize', 11);

    y_est = polyval(f,x);
    hold on
    plot(x,y_est,'k--','LineWidth',2)
    hold off
    
    xlabel(Xdata{i})
    if i==1 || i==3
        ylabel('Correlation to Behaviour')';
    end
    [r,p]=corr(data_Decod(:,i),data_Corr);
    text(max(data_Decod(:,i)),-0.2,['r = ',num2str(r,'%.2f')],'FontSize',10)
    text(max(data_Decod(:,i)),-0.25,['p = ',num2str(p,'%.2f')],'FontSize',10)
end
