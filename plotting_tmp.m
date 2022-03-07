clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
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
    plot(smooth(nanmean(accuracies(Feature,:,:),3),20),'linewidth',3)
    hold on;
end
legend 'mean' 'median' 'variance' 'skewness' 'Kurtosis' 'LZ' 'Higuchi'

% subplot(3,1,2)
figure
for Feature=[9 11:13 18:20]
    plot(smooth(nanmean(accuracies(Feature,:,:),3),20),'linewidth',3)
    hold on;
end
legend 'Katz' 'Hurt' 'SamEntropy' 'ApprEntropy' 'Autocor' 'HjorthComp' 'HjorthMob'

% subplot(3,1,3)
figure;
for Feature=[21:27]
    plot(smooth(nanmean(accuracies(Feature,:,:),3),20),'linewidth',3)
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
for Subject=[1:10]
    %% all channels, 4-class, variance explained: above chance (50.5) decoding: 52.7
    %     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat.mat']);
    
    %% all channels, 4-class, variance explained PCA: running, no bias observed till subj 3
    %     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_PCA.mat']);
    
    
    
    %% all channels, 2-class, sequential: almost no bias, 53.5%
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_Seq.mat']);
    
    %% all channels, 2-class, sequential PCA: running: no bias and 54% very good but better if higher
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_Seq_PCA.mat']);
    
    %% all channels, 2-class, sequential PCA always Mean included: very much like mean alone and better than the above case
%            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_seq_PCA_MeanIs.mat']);
    
    %% all channels, 2-class, sequential PCA always Mean included 2: using keepin argument in the function: no difference almost
%            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_seq_PCA_MeanIs2.mat'],'accuracy','Sel_feat');
    %% all channels, 2-class, sequential PCA Unzero: not good
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_seq_PCA_Unz.mat'],'accuracy','Sel_feat');
    



    %% Every channels, 2-class, correlation: huge bias up to 55% and reaches 57%
    %     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_SelFeat_Cor.mat']);
    
    %% all channels, 2-class, variance explained: bias to 52% and reaches 54.5%
    %     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_Varnc.mat']);
    
    %% all channels, 2-class, variance explained PCA
    %     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_Varnc_PCA.mat']);
    
    %% all channels, 2-class, variance explained PCA UnZ socred: no significant difference from above
%         load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_Varnc_PCA_UnZ.mat']);
    
    
     %% all channels, 2-class, NCA (5-feats): very beautiful, no bias, but 53
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_NCA.mat'],'accuracy','Sel_feat');
    
     %% all channels, 2-class, NCA (5-feats) PCA: like above but 52.5
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_NCA_PCA.mat'],'accuracy','Sel_feat');
    

    %% all channels, 2-class, Relief (5-feats): no bias, good but not very strong, 52.5
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_relieff.mat'],'accuracy','Sel_feat');
    
     %% all channels, 2-class, Relief (5-feats) PCA: like above even weaker 52.2
%     load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_relieff_PCA.mat'],'accuracy','Sel_feat');
    
    accuracies(:,Subject)=nanmean(accuracy,2);
    try
        Feat_Select_all(:,:,Subject)=squeeze(nanmean(Sel_feat,2));
    catch
        Feat_Select_all(:,:,Subject)=squeeze((Sel_feat));
    end
end

% subplot(3,1,1)
figure;
plot(smooth(nanmean(accuracies(:,:),2),20),'linewidth',3)

figure;
try
    imagesc(squeeze(nanmean(Feat_Select_all,3)))
catch
    imagesc(squeeze(nanmean(Feat_Select_all,3)))
end