clc;
clear all;
% close all;
Dataset=1;

Subjects=[1:10];
bands=[1]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:231];
Fs=1000;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
%% These features provide 1 number for each trial and each channel
for band=bands   % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
    
    if band==2
        lowband=0.5;
        highband=4;
    elseif band==3
        lowband=4;
        highband=8;
    elseif band==4
        lowband=8;
        highband=12;
    elseif band==5
        lowband=12;
        highband=16;
    elseif band==6
        lowband=25;
        highband=200;
    end
    
    for Subject=Subjects
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_mrmr_PCA_5ms.mat']);
        accuracies(Subject,:)=squeeze(nanmean(accuracy,2));
        Sel_feats(:,:,Subject)=squeeze(nanmean(Sel_feat,2));
    end
end
plot(nanmean(accuracies))
figure
imagesc(squeeze(nanmean(Sel_feats,3)))

% plot(squeeze(nanmean(Sel_feats,3)))