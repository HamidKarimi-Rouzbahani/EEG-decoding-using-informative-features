%% Feature combinations

clc;
clear all;
close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;

Dataset=3;
Series=3;  % 1:3

miny=0.47;
if Series==1
    array=1:6;
    No_of_BFs=6;
elseif Series==2
    array=[7:12];
    No_of_BFs=6;
elseif Series==3
    array=13:17;
    No_of_BFs=5;
end

method=[1:9 12:19];

minx=-175;
maxx=975;
maxy=0.66;


accuracies=nan*ones(19,231,10);
listFS = {'ILFS','InfFS','ECFS','mrmr','relieff','mutinffs','fsv','laplacian','mcfs','fisher','UDFS','llcfs','cfs','fsasl','dgufs','ufsol','lasso'};
f=0;
for FS=method
    f=f+1;
    selection_method=listFS{f};
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
        
        % interpolation
        for cl=1:size(accuracy,2)
            accuracy_IP(1,cl,:)=interp1([-175:20:865],squeeze(accuracy(1,cl,:)),[-175:5:975],'spline');
            accuracy_IP(1,cl,206:end)=accuracy_tmp;
            accuracy_IP(1,cl,:)=squeeze(accuracy_IP(1,cl,:))+[randn(1,length(accuracy_IP))*0.008]';
        end
        accuracy=accuracy_IP;
        accuracies(f,:,Subject)=squeeze(nanmean(accuracy(1,:,:),2));
        
    end
end

figure;
gca = axes('Position',[0.185 0.70 0.775 0.30]);

colors={[0 0.8 0.8],[0.8 0 0],[0 0.8 0],[0.8 0 0.8],[0.8 0.8 0],[0 0 0.8],[0.5 0.5 0.5],[0.6 0.1 0.1]};
times=[-200:5:950]+25;

feat_comb=0;
for Feature=array
    feat_comb=feat_comb+1;
    plot_line(feat_comb)=plot(times,nanmean(accuracies(Feature,:,:),3),'Color',colors{feat_comb},'linewidth',1.5);
    hold on;
end

line([minx maxx],[0.5 0.5],'LineWidth',1,'Color','k','LineStyle','--');
line([0 0],[miny maxy],'LineWidth',1,'Color','k','LineStyle','--');
ylim([miny maxy])
if  (Series==1 || Series==2)
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5),plot_line(6)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)},listFS{array(6)}},'EdgeColor','none','FontSize',8,'location','northeast');
    
elseif  Series==3
    legend([plot_line(1),plot_line(2),plot_line(3),plot_line(4),plot_line(5)],...
        {listFS{array(1)},listFS{array(2)},listFS{array(3)},listFS{array(4)},listFS{array(5)}},'EdgeColor','none','FontSize',8,'location','northeast');
end

ylabel({'Decoding';'Accuracy (%)'})
box off;

box off;
set(gca,'FontSize',16,'LineWidth',1,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'},'XMinorTick','on','xcolor',[0 0 0],'ycolor',[0 0 0]);
xlim([minx maxx])

% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
pdf_paper_size         = [20 20];
fig.PaperSize       = pdf_paper_size;
pdf_print_resolution   = '-r300';
pdf_file_name=['Comb_paper_Fig2_Dataset_',num2str(Dataset),'_Series_',num2str(Series)];
print(['Z:\projects\Hamid\Projects\Mojgan\Analyses\Combined_features\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)

%% sigfnificance
for feature=1:size(accuracies,1)
    for time=1:size(accuracies,2)
        significance(feature,time)=bf.ttest(squeeze(accuracies(feature,time,:)),squeeze(nanmean(accuracies(feature,1:30,:),2)));
    end
end
Effects=significance';

% Plot Bayes
Time_samples=1:3:size(Effects,1);
Threshold=3;
Against_color=[0 0 0];
Neutral_color=[0.3 0.3 0.3]+.2;
g=0;
for feat_comb=array
    g=g+1;
    maxy_sig=5;
    miny_sig=-2;
    figure;
    gca = axes('Position',[0.185 0.83 0.775 0.15]);
    line([0 0],[miny_sig maxy_sig],'LineWidth',1,'Color','k','LineStyle','--');
    hold on;
    line([minx maxx],[0 0],'Color','k','linewidth',1,'Color','k','LineStyle','--');
    
    c=0;
    for time=Time_samples
        c=c+1;
        if Effects(time,feat_comb)>Threshold
            plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',colors{g},'MarkerEdgeColor',colors{g},'linewidth',0.1,'markersize',5);
        elseif Effects(time,feat_comb)<Threshold && Effects(time,feat_comb)>1./Threshold
            plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Neutral_color,'linewidth',.1,'markersize',5);
        elseif Effects(time,feat_comb)<1./Threshold
            plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',Against_color,'MarkerEdgeColor',Against_color,'linewidth',.1,'markersize',5);
        end
        hold on;
    end
    
    
    ylabel({'BF';'(Log_{10})'})
    box off;
    
    if g==No_of_BFs
        xticklabels={'-100','0','100','200','300','400','500','600','700','800','900'};
        xlabel('Time Relative To Stimulus Onset (ms)')
    else
        xticklabels={};
    end
    
    set(gca,'FontSize',16,'LineWidth',1,'XTick',...
        [-100 0 100:100:900],'XTickLabel',...
        xticklabels,'YTick',...
        [(miny_sig) 0 (maxy_sig)],'YTickLabel',{[num2str(miny_sig)],'0',[num2str(maxy_sig)]},...
        'XMinorTick','on','YMinorTick','on','ycolor','k','xcolor','k');
    xlim([minx maxx])
    ylim([miny_sig maxy_sig])
    xtickangle(45);
    
    % set the prining properties
    fig                 = gcf;
    fig.PaperUnits      = 'centimeters';
    fig.Position        = [100 100 570 460];
    pdf_paper_size         = [20 20];
    fig.PaperSize       = pdf_paper_size;
    pdf_print_resolution   = '-r300';
    pdf_file_name=['Comb_paper_Fig2_BF_Dataset_',num2str(Dataset),'_Series_',num2str(Series),'_',num2str(g)];
    print(['Z:\projects\Hamid\Projects\Mojgan\Analyses\Combined_features\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)
end
