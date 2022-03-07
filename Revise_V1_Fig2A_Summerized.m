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
%     'Mean Freq','Med freq','Avg Freq','SEF 95%','Pw MdFrq','Phs MdFrq','Signal Pw','Cross Corr','Wavelet',...
%     'Hilb Amp','Hilb Phs','CNN','Orig Mag','Amp Lock','Phs Lock'};

Features={'Mean','P1','N1','P2a','P2b','Wavelet','Orig Mag'};
chosen_features_order=[2 14:17 28 34];

band=1;
windoww=1;
Dataset=3;
xtic=[1:1.8:13];

MarkSize=4;

colors={[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8]};
significance=nan*ones(length(bands),length(chosen_features_order));


figure;
gca = axes('Position',[0.16 0.35 0.83 0.1]); %x position, y position, x size and y size

for band=bands
    for Subject=1:10
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
    end
    for feature=1:length(chosen_features_order)
        significance(band,feature)=bf.ttest2(accuracies(feature,:),randn(100,1)*std(accuracies(feature,:))+0.5);
    end
end
% Bayes stats againts chance
for feature=1:length(chosen_features_order)
    Effects=significance';
    for e=1:size(Effects,2)
        if Effects(feature,e)>10
            Bayes(feature,e)=2.5;
        elseif Effects(feature,e)>3 && Effects(feature,e)<=10
            Bayes(feature,e)=1.5;
        elseif Effects(feature,e)>1 && Effects(feature,e)<=3
            Bayes(feature,e)=0.5;
        elseif Effects(feature,e)<1 && Effects(feature,e)>=1/3
            Bayes(feature,e)=-0.5;
        elseif Effects(feature,e)<1/3 && Effects(feature,e)>=1/10
            Bayes(feature,e)=-1.5;
        elseif Effects(feature,e)<1/10
            Bayes(feature,e)=-2.5;
        end
    end
end
for feature=1:length(xtic)
    line([xtic(feature) xtic(feature)]+0.88,[0.39 0.43],'Color','k','linestyle',':','linewidth',1);
    hold on;
end

Baseline=0.45;
steps=0.005;
distans=1.5; % times step
for band=1:size(Bayes,2)
    for feature=1:size(Bayes,1)
        hold on;
        if Bayes(feature,band)==-0.5 || Bayes(feature,band)==0.5
            plots(band)=plot(xtic(feature)+band/4,Bayes(feature,band).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{band},'linewidth',2,'markersize',MarkSize);
        elseif Bayes(feature,band)~=0
            plots(band)=plot(xtic(feature)+band/4,Bayes(feature,band).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{band},'Color',colors{band},'linewidth',2,'markersize',MarkSize);
        end
    end
    baseline_temp=Baseline-(3*2+distans)*steps;
    line([-5 max(xtic)+5],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',0.5);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',0.5);
end

set(gca,'FontSize',16,'DefaultTextFontName','Calibri','XTick',...
    [xtic]+0.88,'XTickLabel',Features,'YTick',...
    [1],'YTickLabel',{''});


ylabel({'Bayes';'Factors'},'FontAngle','italic')
xtickangle(45);
xlim([0.6 14])
ylim([0.39 0.43])
box off;
%% Bar plots
gca = axes('Position',[0.16 0.45 0.83 0.55]);%x position, y position, width and height

for band=bands
    for Subject=1:10
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
    end
    Bars{band}=bar(xtic+band/4,nanmean(accuracies,2),0.14,'facecolor',colors{band},'edgecolor','none');
    hold on;
    errorbar(xtic+band/4,nanmean(accuracies,2),nanstd(accuracies')./sqrt(10),'linewidth',2,'color',colors{band},'CapSize',0,'LineStyle','none')
end
line([-5 max(xtic)+5],[0.5 0.5],'color','k','linestyle','--')
set(gca,'FontSize',16,'DefaultTextFontName','Calibri','XTick',[xtic]+0.88,'XTickLabel',{''},'YTick',[0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'});
box off;
ylabel({'Decoding';'Accuracy (%)'})
xlim([0.6 14])
ylim([0.45 0.68]);
% legend([Bars{1},Bars{2},Bars{3},Bars{4},Bars{5},Bars{6}],{'Broad','Delta(\delta)','Theta(\theta)','Alpha(\alpha)',...
%     'Beta(\beta)','Gamma(\gamma)'},'location','northwest','EdgeColor',[1 1 1],'FontSize',11);
% 
legend([Bars{1},Bars{2},Bars{3},Bars{4},Bars{5},Bars{6}],{'Broad','Delta','Theta','Alpha',...
    'Beta','Gamma'},'location','northwest','EdgeColor',[1 1 1],'FontSize',11);


%% Cross significance 1
figure;
Features={'Mean','P1','N1','P2a','P2b','Wavelet','Orig Mag'};
chosen_features_order=[2 14:17 28 34];

Significanc_mat=nan(6,6,7);
for band1=bands
    for band2=bands
        for Subject=1:10
            load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band1},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
            if band1>1
                accuracy(21:26,:,:)=nan;
            end
            accuracies1(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
            load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band2},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
            if band2>1
                accuracy(21:26,:,:)=nan;
            end
            accuracies2(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
        end
        for feature=1:7
            Significanc_mat(band1,band2,feature)=bf.ttest2(accuracies1(feature,:),accuracies2(feature,:));
            if band1==band2
                Significanc_mat(band1,band2,feature)=0.01;
            end
            if band1<band2
                Significanc_mat(band1,band2,feature)=1;
            end
        end
    end
end

f=0;
for feature=1:length(chosen_features_order)
    for band1=bands
        for band2=bands
            if Significanc_mat(band1,band2,feature)>10
                Bayes_mat(band1,band2,feature)=6;
            elseif Significanc_mat(band1,band2,feature)>3 && Significanc_mat(band1,band2,feature)<=10
                Bayes_mat(band1,band2,feature)=5;
            elseif Significanc_mat(band1,band2,feature)>1 && Significanc_mat(band1,band2,feature)<=3
                Bayes_mat(band1,band2,feature)=4;
            elseif Significanc_mat(band1,band2,feature)<1 && Significanc_mat(band1,band2,feature)>=1/3
                Bayes_mat(band1,band2,feature)=3;
            elseif Significanc_mat(band1,band2,feature)<1/3 && Significanc_mat(band1,band2,feature)>=1/10
                Bayes_mat(band1,band2,feature)=2;
            elseif Significanc_mat(band1,band2,feature)<1/10 && Significanc_mat(band1,band2,feature)~=0
                Bayes_mat(band1,band2,feature)=1;
            elseif Significanc_mat(band1,band2,feature)==0
                Bayes_mat(band1,band2,feature)=0;
            end
        end
    end
    
    f=f+1;
    %     if feature==6 || feature==19 || feature==20
    %         f=f+1;
    %     end
    subplot_tmp=subplot(1,7,f);
    hold(subplot_tmp,'on');
    image(Bayes_mat(:,:,feature),'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    set(subplot_tmp,'CLim',[0 6],'DataAspectRatio',[1 1 1],'FontSize',12,'FontName','Calibri');
    
    xticks(1:length(bands))
    yticks(1:length(bands))
%     xticklabels({'Brd','\delta','\theta','\alpha','\beta','\gamma'})
%     yticklabels({'Brd','\delta','\theta','\alpha','\beta','\gamma'})
    xticklabels([])
    yticklabels([])
    ytickangle(45)
    xtickangle(45)
    title(Features{1,feature})
    
    colormap_mine=parula(6);
    colormap_mine=vertcat([1 1 1],colormap_mine);
    colormap(colormap_mine);
end
%% cross-condition 2
figure;
for band=1:length(bands)
    clearvars Significanc_mat Bayes_mat
    chosen_features_order=[2 14:17 28 34];
    Features={'Mean','P1','N1','P2a','P2b','Wavelet','Orig Mag'};

    accuracies=nan*ones(length(chosen_features_order),10);
    
    for Subject=1:10
        Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
        load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
    end
    
    Significanc_mat=nan(length(chosen_features_order));
    for feature1=1:length(chosen_features_order)
        for feature2=1:length(chosen_features_order)
            Significanc_mat(feature1,feature2)=bf.ttest2(accuracies(feature1,:),accuracies(feature2,:));
            
            if feature1==feature2
                Significanc_mat(feature1,feature2)=0.01;
            end
            if feature1<feature2
                Significanc_mat(feature1,feature2)=0;
            end
            
            if Significanc_mat(feature1,feature2)>10
                Bayes_mat(feature1,feature2)=6;
            elseif Significanc_mat(feature1,feature2)>3 && Significanc_mat(feature1,feature2)<=10
                Bayes_mat(feature1,feature2)=5;
            elseif Significanc_mat(feature1,feature2)>1 && Significanc_mat(feature1,feature2)<=3
                Bayes_mat(feature1,feature2)=4;
            elseif Significanc_mat(feature1,feature2)<1 && Significanc_mat(feature1,feature2)>=1/3
                Bayes_mat(feature1,feature2)=3;
            elseif Significanc_mat(feature1,feature2)<1/3 && Significanc_mat(feature1,feature2)>=1/10
                Bayes_mat(feature1,feature2)=2;
            elseif Significanc_mat(feature1,feature2)<1/10 && Significanc_mat(feature1,feature2)~=0
                Bayes_mat(feature1,feature2)=1;
            elseif Significanc_mat(feature1,feature2)==0
                Bayes_mat(feature1,feature2)=0;
            end
            
            
        end
    end
    subplot_tmp=subplot(1,6,band);
    hold(subplot_tmp,'on');
    image(Bayes_mat,'Parent',subplot_tmp,'CDataMapping','scaled');
    axis(subplot_tmp,'tight');
    axis(subplot_tmp,'ij');
    if band==1
        set(subplot_tmp,'CLim',[0 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    else
        set(subplot_tmp,'CLim',[0 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
    end
    xticks(1:length(chosen_features_order))
    yticks(1:length(chosen_features_order))
%     xticklabels(Features)
%     yticklabels(Features)
    xticklabels([])
    yticklabels([])
    Bands={'Broad','Delta','Theta','Alpha','Beta','Gamma'};
    title(Bands{band},'FontSize',14)
    ytickangle(45)
    xtickangle(45)
    
    colormap_mine=parula(6);
    colormap_mine=vertcat([1 1 1],colormap_mine);
    colormap(colormap_mine);
end
% c = colorbar ('FontSize',16,'Limits',[0.9 6],'YTick',[1:6],'YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
% c.Ticks = [1+0.4 2.1 3 3.85 4.7 5.55]; %Create ticks
% c.Label.String = 'Bayes Factors';
