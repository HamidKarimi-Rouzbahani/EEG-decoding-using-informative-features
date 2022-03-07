clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
windoww=1;

Dataset=3;
% for Dataset=1:3
%    subplot(1,3,Dataset) 
    
    Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
        'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
        'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq',...
        'Cros Cor','Wavelet','Hilb Amp','Hilb Phs','Samples'};
    features=[2:8 9 11:13 18 19 20 27 21:26 32 28:30 34];
    feats=1;
    
    if feats==1
        chosen_features=features(2); %Mean
    else
        chosen_features=features(23); %Wavelet
    end
    minx=-175;
    maxx=975;
    maxy=0.68;
    miny=0.45;
    
    band=1;
    accuracies_broad=nan*ones(1,231,10);
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_cmplt.mat'],'accuracy');
        accuracies_broad(:,:,Subject)=nanmean(accuracy(chosen_features,:,:),2);
    end
    
    band=3;
    accuracies_theta=nan*ones(1,231,10);
    for Subject=[1:10]
        if feats==1
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'.mat'],'accuracy');
        else
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_Wave.mat'],'accuracy');
        end
        accuracies_theta(:,:,Subject)=nanmean(accuracy(chosen_features,:,:),2);
    end
    
    accuracies(1,:,:)=squeeze(accuracies_broad);
    accuracies(2,:,:)=squeeze(accuracies_theta);
    %% Plotting
    
%    subplot(1,3,Dataset) 
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
    
    
    legend([plot_line{1, 1}.mainLine,plot_line{1, 2}.mainLine],{'Broad band','Theta band'},'EdgeColor','w','FontSize',20);
    
    
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
    if feats==1
        title('Mean')
    else
        title('Wavelet')
    end
    
% end