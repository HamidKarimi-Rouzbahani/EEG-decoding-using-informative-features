clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;

Dataset=1;
series=5;  % 1:4


if series==1
    array=1:5;
    miny=0.28;
elseif series==2
    array=6:14;
    miny=0.12;
elseif series==3
    array=15:21;
    miny=0.2;
elseif series==4
    array=22:26;
    miny=0.28;
elseif series==5
    array=27:28;
    miny=0.28;    
end
    
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 34 36 37];    
chosen_features=features(array); %


Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
'Cros Cor','Wavelet','Hilb Amp','Hilb Phs','Samples'};
minx=-175;
maxx=975;
maxy=0.66;

accuracies=nan*ones(37,53,10);
for Subject=[1:10]
    load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
    accuracies(:,:,Subject)=nanmean(accuracy,2);
end

%% Plotting

figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:20:850]+25;

p=0;
for Feature=chosen_features
    p=p+1;
    plot_line(p)=plot(times,smooth(nanmean(accuracies(Feature,:,:),3),1),'Color',colors{p},'linewidth',3);
    hold on;
end

line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');
