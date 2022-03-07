clear all;
close all;
bands=[1]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:8];
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=3;

features=[2 3];
accuracies=nan*ones(length(features),231,10);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies(:,:,Subject)=squeeze(nanmean(accuracy(features,:,:),2));
end

plot(squeeze(nanmean(accuracies,3))')

figure;
features=[28:30 32 34:35];
accuracies=nan*ones(length(features),231,10);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'mult.mat'],'accuracy');
    accuracies(:,:,Subject)=squeeze(nanmean(accuracy(features,:,:),2));
end

plot(squeeze(nanmean(accuracies,3))')
%% Completing the datasets
clear all;
close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=3;

% features=[2 3];
accuracy_tmp=nan*ones(35,15,231,2);
for Subject=[1:10]
    clearvars -except windoww Dataset band Bands Subject accuracy_tmp Datasets

    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracy_tmp(:,:,:,1)=accuracy;    
    
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'mult.mat'],'accuracy');
    accuracy_tmp(:,:,:,2)=accuracy;    
    accuracy_final=nansum(accuracy_tmp,4);
    accuracy_final(accuracy_final==0)=nan;
    accuracy=accuracy_final;
    save(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_cmplt.mat'],'accuracy');
end

plot(squeeze(nanmean(accuracy_final([2 28:35],:,:),2))')
%%

clear all;
close all;
bands=[1]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:8];
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=1;

feature=[1];
accuracies=nan*ones(53,10);
Sel_feats=nan*ones(53,26,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','rfe','L0','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
for Subject=[1:9]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',listFS{feature},'_PCA_cmplt.mat'],'accuracy','Sel_feat');
    accuracies(:,Subject)=squeeze(nanmean(accuracy,2));
    Sel_feats(:,:,Subject)=squeeze(nanmean(Sel_feat,2));
end

plot(nanmean(accuracies,2))
figure;
imagesc(nanmean(Sel_feats,3))
