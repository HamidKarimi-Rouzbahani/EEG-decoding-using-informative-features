%% Feature combinations

clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=2;

series=1;  % 1:3
if series==1
    array=1:17;
    miny=0.24;
elseif series==2
    array=[7:12];
    miny=0.24;
elseif series==3
    array=13:17;
    miny=0.28;
end
    
method=[1:9 12:19];    

minx=-175;
maxx=975;
maxy=0.64;


accuracies=nan*ones(19,53,10);
Feat_Select_all=nan*ones(53,26,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};

feature_maps=nan*ones(17,53,26,10);
f=0;
for feat=method
    f=f+1;
    selection_method=listFS{f};
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy','Sel_feat');
        accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));   
        Features_used(:,:,f,Subject)=squeeze(nanmean(Sel_feat,2));
    end
end

figure;
% colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:20:950]+25;

p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(smooth(nanmean(accuracies(Feature,:,:),3),3),'linewidth',3);
    hold on;
end
figure;
for i=1:17
    subplot(6,3,i)
    imagesc(squeeze(nanmean(Features_used(:,:,i,:),4))');
end
% line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
% line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');