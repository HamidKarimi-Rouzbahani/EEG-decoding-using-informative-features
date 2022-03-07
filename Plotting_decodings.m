clc;
clear all;
close all;



% for Subject=1:10
%     load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'.mat'],'accuracy','True_Predicted_labels');
%     accuracies(:,:,:,Subject)=accuracy;
%     
% end
% 
% plot(nanmean(nanmean(nanmean(accuracies,4),3),2))
% hold on;
% for Subject=1:4
%     load(['Decoding_Accuracy_DS_Mine_LDA_Subject_',num2str(Subject),'.mat'],'accuracy','True_Predicted_labels');
%     accuracies(:,:,:,Subject)=accuracy;
%     
% end
% plot(nanmean(nanmean(nanmean(accuracies(:,:,:,1:4),4),3),2))
% 

for Subject=1:10
%     Decoding_Accuracy_ALL_Subject_10VahabsDS_filtered_Beta
    load(['Decoding_Accuracy_ALL_Subject_',num2str(Subject),'VahabsDS_filtered_Gamma.mat'],'accuracy','True_Predicted_labels');
    accuracies(:,:,:,Subject)=accuracy;
    
end

plot(nanmean(nanmean(nanmean(accuracies,4),3),2))
hold on;
for Subject=1:4
    load(['Decoding_Accuracy_DS_Vahab_LDA_Subject_',num2str(Subject),'.mat'],'accuracy','True_Predicted_labels');
    accuracies(:,:,:,Subject)=accuracy;
    
end
plot(nanmean(nanmean(nanmean(accuracies(:,:,:,1:4),4),3),2))

