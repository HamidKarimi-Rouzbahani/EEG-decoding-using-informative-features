%% plotting
clc;
clear all;
% close all;
bands=[1:6]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[2:8];
Datasets={'Mine','Vhab','Stfd'};
% DS1: Animal Car Face Plane
% DS2: Object Face Animal Fruit
% DS3: AnimalBody AnimalFace ruitVegetable HumanBody HumanFace InanimateObject

addpath('K:\MQ_Analysis_PC_13_3_2020\Hamid\Mojgan_analyses\Analyses\bayesFactor-master');
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
Features={'Baseline','Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx',...
    'Higuchi FD','Katz FD','Hurst Exp','Sample Ent','Apprx Ent',...
    'P1','N1','P2a','P2b','Autocorr','Hjorth Cmp','Hjorth Mob','Signal Pw',...
    'Mean Freq','Med freq','Avg Freq','SEF 95%','Pw MdFrq','Phs MdFrq','Wavelet',...
    'Hilb Amp','Hilb Phs','Cros Cor','CNN','Samples','AutoCorM'};
chosen_features_order=[1:9 11:20 27 21:26 28:30 32:35];

band=1;
windoww=1;
Dataset=3;
xtic=[1:1.8:1000];
xtic=xtic(1:33);
xtic(21:26)=xtic(21:26)+1;
xtic(27:33)=xtic(27:33)+2;

colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5]};
significance=nan*ones(length(bands),length(chosen_features_order));


figure;
gca = axes('Position',[0.13 0.131 0.775 0.2]);

for windoww=Windows
    for Subject=1:10
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
    end
    for feature=1:length(chosen_features_order)
        significance(windoww,feature)=bf.ttest2(accuracies(feature,:),randn(100,1)*std(accuracies(feature,:))+0.5);
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
    line([xtic(feature) xtic(feature)]+1.08,[0.39 0.43],'Color','k','linestyle',':','linewidth',1);
    hold on;
end

Baseline=0.45;
steps=0.005;
distans=1.5; % times step
for windoww=1:size(Bayes,2)
    for feature=1:size(Bayes,1)
        hold on;
        if Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==0.5
            plots(windoww)=plot(xtic(feature)+windoww*(7/8)/4,Bayes(feature,windoww).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{windoww},'linewidth',2,'markersize',10);
        elseif Bayes(feature,windoww)~=0
            plots(windoww)=plot(xtic(feature)+windoww*(7/8)/4,Bayes(feature,windoww).*steps+Baseline-(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{windoww},'Color',colors{windoww},'linewidth',2,'markersize',10);
        end
    end
    baseline_temp=Baseline-(3*2+distans)*steps;
    line([-5 max(xtic)+5],[baseline_temp baseline_temp],'linestyle','--','Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-2*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]-3*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+2*steps,'Color','k','linewidth',1);
    line([-5 max(xtic)+5],[baseline_temp baseline_temp]+3*steps,'Color','k','linewidth',1);
end

set(gca,'FontSize',18,'FontName','Calibri','XTick',...
    [xtic]+1.08,'XTickLabel',Features,'YTick',...
    [1],'YTickLabel',{''});


ylabel({'Bayes';'Factors'})
xtickangle(45);
xlim([2.6 63])
ylim([0.39 0.43])
box off;
%% Bar plots
gca = axes('Position',[0.13 0.05 0.775 0.794]);
for windoww=Windows
    for Subject=1:10
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
    end
    Bars{windoww}=bar(xtic+windoww*(7/8)/4,nanmean(accuracies,2),0.14*7/8,'facecolor',colors{windoww},'edgecolor','none');
    hold on;
    errorbar(xtic+windoww*(7/8)/4,nanmean(accuracies,2),nanstd(accuracies')./sqrt(10),'linewidth',2,'color',colors{windoww},'CapSize',0,'LineStyle','none')
end
line([-5 max(xtic)+5],[0.5 0.5],'color','k','linestyle','--')
set(gca,'FontSize',18,'FontName','Calibri','XTick',[xtic]+1.08,'XTickLabel',{''},'YTick',[0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'});
box off;
ylabel('Decoding Accuracy (%)')
xlim([2.6 63])
ylim([0.45 0.68]);
legend([Bars{2},Bars{3},Bars{4},Bars{5},Bars{6},Bars{7},Bars{8}],{'-200:0 ms','0:100 ms',...
    '100:200 ms','200:300 ms','300:400 ms','400:500 ms','500:600 ms'},'location','northwest','EdgeColor',[1 1 1]);
ccc
%% Cross significance
Windows=1:7;
figure;
for windoww1=Windows
    for windoww2=Windows
        for Subject=1:10
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww1+1),'_Subject_',num2str(Subject),'.mat'],'accuracy');
%             if band1>1
%                 accuracy(21:26,:,:)=nan;
%             end
            accuracies1(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
            
            load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww2+1),'_Subject_',num2str(Subject),'.mat'],'accuracy');
%             if band2>1
%                 accuracy(21:26,:,:)=nan;
%             end
            accuracies2(:,Subject)=nanmean(nanmean(accuracy(chosen_features_order,:,:),2),3);
        end
        for feature=1:33
            Significanc_mat(windoww1,windoww2,feature)=bf.ttest2(accuracies1(feature,:),accuracies2(feature,:));
            if windoww1==windoww2
                Significanc_mat(windoww1,windoww2,feature)=1;
            end
        end
    end
end
f=0;
for feature=1:length(chosen_features_order)
    for windoww1=Windows
        for windoww2=Windows
            if Significanc_mat(windoww1,windoww2,feature)>10
                Bayes_mat(windoww1,windoww2,feature)=6;
            elseif Significanc_mat(windoww1,windoww2,feature)>3 && Significanc_mat(windoww1,windoww2,feature)<=10
                Bayes_mat(windoww1,windoww2,feature)=5;
            elseif Significanc_mat(windoww1,windoww2,feature)>1 && Significanc_mat(windoww1,windoww2,feature)<=3
                Bayes_mat(windoww1,windoww2,feature)=4;
            elseif Significanc_mat(windoww1,windoww2,feature)<1 && Significanc_mat(windoww1,windoww2,feature)>=1/3
                Bayes_mat(windoww1,windoww2,feature)=3;
            elseif Significanc_mat(windoww1,windoww2,feature)<1/3 && Significanc_mat(windoww1,windoww2,feature)>=1/10
                Bayes_mat(windoww1,windoww2,feature)=2;
            elseif Significanc_mat(windoww1,windoww2,feature)<1/10
                Bayes_mat(windoww1,windoww2,feature)=1;
            end
        end
    end
    
    if feature>1 && (feature<13 || feature>16) && feature~=31
        f=f+1;
        if feature<33
            subplot_tmp=subplot(3,9,f);
            hold(subplot_tmp,'on');
            image(Bayes_mat(:,:,feature),'Parent',subplot_tmp,'CDataMapping','scaled');
            axis(subplot_tmp,'tight');
            axis(subplot_tmp,'ij');
            set(subplot_tmp,'CLim',[1 6],'DataAspectRatio',[1 1 1],'FontSize',10,'FontName','Calibri');
            
            xticks(1:length(Windows))
            yticks(1:length(Windows))
            xticklabels({'BL','W1','W2','W3','W4','W5','W6'})
            yticklabels({'BL','W1','W2','W3','W4','W5','W6'})
            ytickangle(45)
            xtickangle(45)
            title(Features{1,feature})
            set(gca,'FontSize',12,'FontName','Calibri');
        else
            figure;
            subplot_tmp=subplot(1,1,1);
            hold(subplot_tmp,'on');
            image(Bayes_mat(:,:,feature),'Parent',subplot_tmp,'CDataMapping','scaled');
            axis(subplot_tmp,'tight');
            axis(subplot_tmp,'ij');
            set(subplot_tmp,'CLim',[1 6],'DataAspectRatio',[1 1 1],'FontSize',12,'FontName','Calibri');
            
            xticks(1:length(Windows))
            yticks(1:length(Windows))
            xticklabels({'BL','W1','W2','W3','W4','W5','W6'})
            yticklabels({'BL','W1','W2','W3','W4','W5','W6'})
            ytickangle(45)
            xtickangle(45)
            colormap(parula(6));
            c = colorbar ('YTickLabel',{'BF<0.1','0.1<BF<0.3','0.3<BF<1','1<BF<3','3<BF<10','BF>10'}) ; %Create Colorbar
            c.Ticks = [1:6] ; %Create ticks
            c.Label.String = 'Bayes Factors';
            title(Features{1,feature})
            set(gca,'FontSize',14,'FontName','Calibri');
        end

    end
end
