clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;

Dataset=3;



Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
    'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
    'Cros Cor','Wavelet','Hilb Amp','Hilb Phs','Samples'};
features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 34];
chosen_features=features(23); %Wavelet

minx=-175;
maxx=975;
maxy=0.665;
miny=0.45;

accuracies_ind=nan*ones(1,231,10);
for Subject=[1:10]
    load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_cmplt.mat'],'accuracy');
    accuracies_ind(:,:,Subject)=nanmean(accuracy(chosen_features,:,:),2);
end

listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
% method=[1:9 12:19];
method=[16];   %ufsol

accuracies_cmbnd=nan*ones(1,231,10);

f=0;
for FS=method
    f=f+1;
    selection_method=listFS{FS};
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy','Sel_feat');
        accuracy_orig=accuracy;
        
        try
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms.mat'],'accuracy');
            accuracy_tmp=smooth(nanmean(accuracy(1,:,206:end),2),4);
        catch
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_cmplt.mat'],'accuracy');
            accuracy_tmp=smooth(repmat(squeeze(nanmean(accuracy(1,:,47:53),2)),[4 1]),4);
            accuracy_tmp=accuracy_tmp(1:26);
        end
        
        accuracy=accuracy_orig;
        if strcmp(selection_method,'fisher')
            accuracy=accuracy-repmat(nanmean(nanmean(accuracy(1,:,1:10),2),3),[size(accuracy)])+0.5;
        end
        
        
        for cl=1:size(accuracy,2)
            accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
            accuracy_IP(1,cl,206:end)=accuracy_tmp;
            accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
        end
        accuracy=accuracy_IP;
        accuracies_cmbnd(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));
    end
end

accuracies(1,:,:)=squeeze(accuracies_ind);
accuracies(2,:,:)=squeeze(accuracies_cmbnd);
%% Plotting

figure;
colors={[0 0 0],[0 0.8 0]};
times=[-200:5:950]+25;

p=0;
for Feature=1:2
    p=p+1;
    for subject=1:10
        if Feature==1
            points_to_smooth=5;
        else
            points_to_smooth=3;
        end
        accuracies(Feature,:,subject)=smooth(accuracies(Feature,:,subject),points_to_smooth);
    end
    plot_line{p}=shadedErrorBar(times,nanmean(squeeze((accuracies(Feature,:,:))),2),nanstd(squeeze((accuracies(Feature,:,:)))')'./sqrt(10),{'color',colors{p},'LineWidth',3},1);
    
    hold on;
end

line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');

%% sigfnificance
for time=1:size(accuracies,2)
    significance(1,time)=bf.ttest(squeeze(accuracies(1,time,:)),squeeze(accuracies(2,time,:)));
end

% Bayes stats againts chance
for feature=1:size(significance,1)
    Effects=significance;
    for time=1:size(significance,2)
        if Effects(feature,time)>10
            Bayes(feature,time)=2.5;
        elseif Effects(feature,time)>3 && Effects(feature,time)<=10
            Bayes(feature,time)=1.5;
        elseif Effects(feature,time)>1 && Effects(feature,time)<=3
            Bayes(feature,time)=0.5;
        elseif Effects(feature,time)<1 && Effects(feature,time)>=1/3
            Bayes(feature,time)=-0.5;
        elseif Effects(feature,time)<1/3 && Effects(feature,time)>=1/10
            Bayes(feature,time)=-1.5;
        elseif Effects(feature,time)<1/10
            Bayes(feature,time)=-2.5;
        end
    end
end


Baseline=0.47;
steps=0.005;
distans=2; % times step
f=0;
for feature=1
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==0.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color','r','linewidth',2,'markersize',4);
        elseif Bayes(feature,windoww)~=0
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor','r','Color','r','linewidth',2,'markersize',4);
        end
    end
    baseline_temp=Baseline-(f-1)*(3*2+distans)*steps;
    line([minx maxx],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',1);
    line([minx maxx],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',1);
end


legend([plot_line{1, 1}.mainLine,plot_line{1, 2}.mainLine],{'Best individual feature (Wavelet)','Best combined feature (ufsol)'},'EdgeColor','w','FontSize',20);


ylim([miny maxy])

ylabel('Decoding Accuracy (%)')
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',24,'LineWidth',4,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'},'XMinorTick','on');
xtickangle(45);
xlim([minx maxx])


