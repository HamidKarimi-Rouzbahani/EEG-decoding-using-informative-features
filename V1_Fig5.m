%% Feature combinations

clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
Dataset=1;

series=1;  % 1:3
if series==1
    array=1:6;
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


accuracies=nan*ones(19,231,10);
Feat_Select_all=nan*ones(231,21,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
features=[2:9 11:13 18:20 27 21:26];    

feature_maps=nan*ones(17,231,21,10);
f=0;
for feat=method
    f=f+1;
    selection_method=listFS{f};
    for Subject=[1:10]
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_sliding_Subject_',num2str(Subject),'_CombFeat_',selection_method,'_PCA_5ms.mat'],'accuracy','Sel_feat');
        accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));   
        feature_maps(f,:,:,Subject)=squeeze(nanmean(Sel_feat,2));
    end
end

figure;
colors={[0 0.8 0.8],[0 0 0],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

p=0;
for Feature=array
    p=p+1;
    plot_line(p)=plot(times,smooth(nanmean(accuracies(Feature,:,:),3),5),'Color',colors{p},'linewidth',3);
    hold on;
end

line([minx maxx],[0.5 0.5],'LineWidth',1.5,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1.5,'Color','k','LineStyle','--');
save('feature_maps.mat','feature_maps')
%% sigfnificance
for feature=1:size(accuracies,1)
    for time=1:size(accuracies,2)
        significance(feature,time)=bf.ttest(squeeze(accuracies(feature,time,:)),squeeze(nanmean(accuracies(feature,1:30,:),2)));
    end
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
for feature=array
    f=f+1;
    hold on;
    for windoww=1:size(Bayes,2)
        if Bayes(feature,windoww)==-0.5 || Bayes(feature,windoww)==0.5
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','Color',colors{f},'linewidth',2,'markersize',4);
        elseif Bayes(feature,windoww)~=0
            plots(f)=plot(times(windoww),Bayes(feature,windoww).*steps+Baseline-(f-1)*(3*2+distans)*steps,'LineStyle','none','marker','o','MarkerFaceColor',colors{f},'Color',colors{f},'linewidth',2,'markersize',4);
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

if series==1 || series==2
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)},listFS{array(6)}},'EdgeColor','w','FontSize',14);

elseif series==3
     legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)}},'EdgeColor','w','FontSize',14);   
end
ylim([miny maxy])
ylabel('Decoding Accuracy')
box off;

xlabel('Time Relative to Stimulus Onset (ms)')
box off;
set(gca,'FontSize',18,'LineWidth',4,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {'-100','0','100','200','300','400','500','600','700','800','900'},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'0.5','0.55','0.6','0.65'},'XMinorTick','on');
xlim([minx maxx])
ccc
%%
clc
clear all;
close all;
load('feature_maps.mat','feature_maps')

Feat_names={'Mean','Median','Variance','Skewness','Kurtosis','LZ Cmplx','Higuchi FD',...
    'Katz FD','Hurst Exp','Sample Ent','Apprx Ent','Autocorr','Hjorth Cmp','Hjorth Mob',...
'Signal Pw','Mean Freq','Med Freq','Avg Freq','SEF 95%','Pw MedFrq','Phs MdFrq'};

listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
times=[-200:5:950]+25;

for i=1:17
   subplot_tmp=subplot(4,5,i);
   hold(subplot_tmp,'on');
   imagesc(squeeze(nanmean(feature_maps(i,:,:,:),4))','Parent',subplot_tmp,'CDataMapping','scaled');
   axis(subplot_tmp,'tight');
   axis(subplot_tmp,'ij');
   set(subplot_tmp,'CLim',[0 1],'FontSize',10,'FontName','Calibri');
   hold on;
   line([36 36],[0.5 21.5],'LineWidth',1.5,'Color','r');
   
   xticks([26 36:20:231])
   yticks([1:21])
   xticklabels({'-100','0','100','200','300','400','500','600','700','800','900'})
   yticklabels(Feat_names)
   ytickangle(45)
   xtickangle(45)
   title(listFS{1,i})
   set(gca,'FontSize',8,'FontName','Calibri');
end

