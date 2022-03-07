clc
clear all
close all

% band=1;
% Subject=1;
Accuracies=nan.*ones(35,10,5);
for band=1:4
    for Subject=[1:3 6:9]
        if band==1
            load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Delta.mat'],'accuracy','True_Predicted_labels');
        elseif band==2
            load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Theta.mat'],'accuracy','True_Predicted_labels');
        elseif band==3
            load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Alpha.mat'],'accuracy','True_Predicted_labels');
        elseif band==4
            load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Beta.mat'],'accuracy','True_Predicted_labels');
        elseif band==5
            load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Gamma.mat'],'accuracy','True_Predicted_labels');
        end
        
        Accuracies(:,Subject,band)=nanmean(nanmean(accuracy,2),3);
        
    end
end
plot(squeeze(nanmean(Accuracies,2)))