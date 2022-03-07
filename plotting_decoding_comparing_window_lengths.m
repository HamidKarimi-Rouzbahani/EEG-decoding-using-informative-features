clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=2;
accuracies1=nan*ones(35,231,10);
accuracies2=nan*ones(35,231,10);
accuracies3=nan*ones(35,221,10);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_5msWind.mat'],'accuracy');
    accuracies1(:,:,Subject)=nanmean(accuracy,2);
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies2(:,:,Subject)=nanmean(accuracy,2);
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_100msWind.mat'],'accuracy');
    accuracies3(:,:,Subject)=nanmean(accuracy,2);
end
figure;
for Feature=[2]
    plot(smooth(nanmean(accuracies1(Feature,:,:),3),4),'linewidth',3)
    hold on;
    plot(smooth(nanmean(accuracies2(Feature,:,:),3),4),'linewidth',3)
    plot(smooth(nanmean(accuracies3(Feature,:,:),3),4),'linewidth',3)
end