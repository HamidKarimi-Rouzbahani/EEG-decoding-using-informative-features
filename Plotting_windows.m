clc;
clear all;
close all;

for Subject=1:10
    Dataset=1;
    windows=[101:200;201:300;301:400;401:500;501:600;601:700;701:800];
    for window=1:7
        if Dataset==1
            load(['Decoding_Accuracy_DS_Mine_Window_',num2str(window),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        elseif Dataset==2
            load(['Decoding_Accuracy_DS_Vahab_Window_',num2str(window),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        elseif Dataset==3
            load(['Decoding_Accuracy_DS_Stanford_Window_',num2str(window),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        end
%         plot(nanmean(nanmean(accuracy,2),3))
%         hold on;
Accuracies(Subject,window,:)=nanmean(nanmean(accuracy,2),3);

    end
    % ccc
end
plot(squeeze(nanmean(Accuracies,1))')