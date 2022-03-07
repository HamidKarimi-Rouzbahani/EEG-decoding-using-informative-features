%% plotting
clc;
clear all;
% close all;
bands=[1:6]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:8];
Datasets={'Mine','Vhab','Stfd'};
% DS1: Animal Car Face Plane
% DS2: Object Face Animal Fruit
% DS3: AnimalBody AnimalFace ruitVegetable HumanBody HumanFace InanimateObject
% chosen_features_order=[1:9 11:20 27 21:26 28:30 32:35];

addpath('K:\MQ_Analysis_PC_13_3_2020\Hamid\Mojgan_analyses\Analyses\bayesFactor-master');
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
% Features={'Baseline','Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx',...
%     'Higuchi FD','Katz FD','Hurst Exp','Sample Ent','Apprx Ent',...
%     'P1','N1','P2a','P2b','Autocorr','Hjorth Cmp','Hjorth Mob',...
%     'Mean Freq','Med freq','Avg Freq','SEF 95%','Pw MdFrq','Phs MdFrq','Signal Pw','Cros Cor','Wavelet',...
%     'Hilb Amp','Hilb Phs','CNN','Samples','AutoCorM'};
% chosen_features_order=[1:9 11:20 21:27 32 28:30 33:35];


Features={'Baseline','Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx',...
    'Higuchi FD','Katz FD','Hurst Exp','Sample Ent','Apprx Ent',...
    'Autocorr','Hjorth Cmp','Hjorth Mob','P1','N1','P2a','P2b',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MdFrq','Phs MdFrq','Cros Cor','Wavelet',...
    'Hilb Amp','Hilb Phs','Samples','AutoCorM'};
chosen_features_order=[1:9 11:13 18:20 14:17 27 21:26 32 28:30 34:37];
% chosen_features_order=[1:9 11 16:18 12:15 19:22 27 23:26 32 28:30 33:35];



band=1;
windoww=1;
Dataset=3;
xtic=[1:1.8:61];
xtic(7:34)=xtic(7:34)+1;
xtic(20:34)=xtic(20:34)+1;
xtic(27:34)=xtic(27:34)+1;

MarkSize=4;

colors={[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8]};
significance=nan*ones(length(bands),length(chosen_features_order));


figure;
gca = axes('Position',[0.13 0.131 0.775 0.2]);

for band=bands
    for Subject=1:10
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracyT=accuracy;
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracyT(36:37,:,:)=accuracy(36:37,:,:);
        
        accuracy=accuracyT;
        save(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
%         if band>1
%             accuracyT(21:26,:,:)=nan;
%         end
%         accuracies(:,Subject)=nanmean(nanmean(accuracyT(chosen_features_order,:,:),2),3);
    end
%     Bars{band}=bar(xtic+band/4,nanmean(accuracies,2),0.14,'facecolor',colors{band},'edgecolor','none');
%     hold on;
%     errorbar(xtic+band/4,nanmean(accuracies,2),nanstd(accuracies')./sqrt(10),'linewidth',2,'color',colors{band},'CapSize',0,'LineStyle','none')
end
ccc
line([-5 max(xtic)+5],[0.5 0.5],'color','k','linestyle','--')
set(gca,'FontSize',24,'DefaultTextFontName','Calibri','XTick',[xtic]+0.88,'XTickLabel',{''},'YTick',[0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'});
box off;
ylabel('Decoding Accuracy (%)')
xlim([2.6 63.5])
ylim([0.45 0.68]);
legend([Bars{1},Bars{2},Bars{3},Bars{4},Bars{5},Bars{6}],{'Broad','Delta(\delta)','Theta(\theta)','Alpha(\alpha)',...
    'Beta(\beta)','Gamma(\gamma)'},'location','northwest','EdgeColor',[1 1 1]);
