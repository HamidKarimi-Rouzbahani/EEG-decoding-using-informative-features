clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=3;
windoww=1;
Dataset=1;
accuracies=nan*ones(35,231,10);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies(:,:,Subject)=nanmean(accuracy,2);
end

% subplot(3,1,1)
figure;
for Feature=[2:8]
    plot(smooth(nanmean(accuracies(Feature,:,:),3),4),'linewidth',3)
    hold on;
end
ccc
legend 'mean' 'median' 'variance' 'skewness' 'Kurtosis' 'LZ' 'Higuchi'

% subplot(3,1,2)
figure
for Feature=[9 11:13 18:20]
    plot(smooth(nanmean(accuracies(Feature,:,:),3),4),'linewidth',3)
    hold on;
end
legend 'Katz' 'Hurt' 'SamEntropy' 'ApprEntropy' 'Autocor' 'HjorthComp' 'HjorthMob'

% subplot(3,1,3)
figure;
for Feature=[21:27]
    plot(smooth(nanmean(accuracies(Feature,:,:),3),4),'linewidth',3)
    hold on;
end
legend 'MeanFreq' 'MedFreq' 'SEF' 'PowMedFreq' 'PhasMedFreq' 'Power'

%% Features combined
clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=1;
accuracies=nan*ones(53,10);
Feat_Select_all=nan*ones(53,21,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','rfe','L0','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
selection_method=listFS{13};

for Subject=[1:10]
    %% all channels, 4-class, variance explained: above chance (50.5) decoding: 52.7
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA.mat'],'accuracy','Sel_feat');
    
    accuracies(:,Subject)=nanmean(accuracy,2);
    try
        Feat_Select_all(:,:,Subject)=squeeze(nanmean(Sel_feat,2));
    catch
        Feat_Select_all(:,:,Subject)=squeeze((Sel_feat));
    end
end

% subplot(3,1,1)
figure;
plot(smooth(nanmean(accuracies(:,:),2),4),'linewidth',3)

figure;
try
    imagesc(squeeze(nanmean(Feat_Select_all,3)))
catch
    imagesc(squeeze(nanmean(Feat_Select_all,3)))
end

%% Features compared
clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=1;
accuracies=nan*ones(53,19,10);
Feat_Select_all=nan*ones(53,21,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','rfe','L0','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
f=0;
for feat=[1]
% for feat=[19:-1:13 9:-1:1]
% for feat=[6]
    selection_method=listFS{feat};
    for Subject=[1:10]
        %% 3 features
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_3Feats.mat'],'accuracy','Sel_feat');
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_3Feats.mat'],'accuracy','Sel_feat');
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_3Feats_v2.mat'],'accuracy','Sel_feat');
        %% 5 features
%                 load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA.mat'],'accuracy','Sel_feat');
%               load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms_2.mat'],'accuracy','Sel_feat');
        %% 5 features non-linear
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_NonLin.mat'],'accuracy','Sel_feat');
        %% 7 features
%                 load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_7Feats.mat'],'accuracy','Sel_feat');
                load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_7Feats.mat'],'accuracy','Sel_feat');
        %% 10 features
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_10Feats.mat'],'accuracy','Sel_feat');

        accuracies(:,feat,Subject)=nanmean(accuracy,2);
        Feat_Select_all(:,:,Subject)=squeeze(nanmean(Sel_feat,2));
    end
end
figure;
plot(nanmean(accuracies,3),'linewidth',3)

% figure;
% imagesc(squeeze(nanmean(Feat_Select_all,3)),[0 1])


%% Features compared
clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=1;
accuracies=nan*ones(231,21,10);
Feat_Select_all=nan*ones(231,21,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','rfe','L0','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
f=0;
feat=[6];
selection_method=listFS{feat};
for Subject=[1:10]
%             load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA.mat'],'accuracy','Sel_feat');        

%% features only 4-6

%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms_2.mat'],'accuracy','Sel_feat');
%% Without mean only for Mutual 
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms_WO_Mean.mat'],'accuracy','Sel_feat');

%% Dataset 3
% load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA.mat'],'accuracy','Sel_feat');

accuracies(:,feat,Subject)=nanmean(accuracy,2);
Feat_Select_all(:,:,Subject)=squeeze(nanmean(Sel_feat,2));

end
% figure;
hold on;
plot(nanmean(accuracies,3),'linewidth',3)

figure;
imagesc(squeeze(nanmean(Feat_Select_all,3)),[0 1])
