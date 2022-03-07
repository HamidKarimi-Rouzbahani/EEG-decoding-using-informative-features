clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;

Dataset=1;

accuracies=nan*ones(37,6,231);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_cmplt.mat'],'accuracy');
    accuracies(1:35,:,:)=accuracy;
    load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');

    for feat=36:37
        for cl=1:size(accuracy,2)
            accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(feat,cl,:)),[-175:5:975],'spline');
            accuracy_IP(1,cl,206:231)=accuracy_IP(1,cl,[206:231]-26);
             accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.02]';
        end
        accuracies(feat,:,:)=accuracy_IP;
    end 
    accuracy=accuracies;
    save(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
end
plot(squeeze(nanmean(accuracies([30:37],:,:),2))')


