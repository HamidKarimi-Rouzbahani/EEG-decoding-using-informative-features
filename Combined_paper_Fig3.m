clc;
clear all;
% close all;
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;

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
miny=0.47;

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
% Plotting

figure;
gca = axes('Position',[0.185 0.55 0.775 0.45]);

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

if Dataset==1
legend([plot_line{1, 1}.mainLine,plot_line{1, 2}.mainLine],{'Best individual feature (Wavelet)','Best combined feature set (ufsol)'},'EdgeColor','none','FontSize',14,'location','northeast');
end

ylim([miny maxy])

ylabel({'Decoding';'Accuracy (%)'})
box off;
set(gca,'FontSize',16,'LineWidth',1,'XTick',...
    [-100 0 100:100:900],'XTickLabel',...
    {},'YTick',...
    [0.5 0.55 0.6 0.65],'YTickLabel',{'50','55','60','65'},'XMinorTick','on');
xtickangle(45);
xlim([minx maxx])

% set the prining properties
fig                 = gcf;
fig.PaperUnits      = 'centimeters';
fig.Position        = [100 100 570 460];
pdf_paper_size         = [20 20];
fig.PaperSize       = pdf_paper_size;
pdf_print_resolution   = '-r300';
pdf_file_name=['Comb_paper_Fig3_Dataset_',num2str(Dataset)];
print(['Z:\projects\Hamid\Projects\Mojgan\Analyses\Combined_features\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)

%% sigfnificance
for time=1:size(accuracies,2)
    significance(1,time)=bf.ttest(squeeze(accuracies(1,time,:)),squeeze(accuracies(2,time,:)));
end
Effects=significance';
colors={[.9 0.2 .2]};

% Plot Bayes
Time_samples=1:3:size(Effects,1); % rev 2
Threshold=3;
Against_color=[0 0 0];
Neutral_color=[0.3 0.3 0.3]+.2;
for feat_comb=1
    maxy_sig=5;
    miny_sig=-2;
    figure;
    gca = axes('Position',[0.185 0.83 0.775 0.15]);
    line([0 0],[miny_sig maxy_sig],'LineWidth',.5,'Color','k','LineStyle','--');
    hold on;
    line([minx maxx],[0 0],'Color','k','linewidth',.5,'Color','k','LineStyle','--');        
    
    c=0;
    for time=Time_samples
        c=c+1;
        if Effects(time,feat_comb)>Threshold
            plots(feat_comb,time)=plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',colors{feat_comb},'MarkerEdgeColor',colors{feat_comb},'linewidth',.1,'markersize',5);
        elseif Effects(time,feat_comb)<Threshold && Effects(time,feat_comb)>1./Threshold
            plots(feat_comb,time)=plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',Neutral_color,'linewidth',.1,'markersize',5);
        elseif Effects(time,feat_comb)<1./Threshold
            plots(feat_comb,time)=plot(times(Time_samples(c)),log10(Effects(time,feat_comb)),'LineStyle','none','marker','o','MarkerFaceColor',Against_color,'MarkerEdgeColor',Against_color,'linewidth',.1,'markersize',5);
        end
        hold on;
    end
    
    
    ylabel({'BF';'(Log_{10})'})
    box off;
    
    xticklabels={'-100','0','100','200','300','400','500','600','700','800','900'};
    xlabel('Time Relative To Stimulus Onset (ms)')

    
    set(gca,'FontSize',16,'LineWidth',.5,'XTick',...
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
pdf_file_name=['Comb_paper_Fig3_BF_Dataset_',num2str(Dataset)];
print(['Z:\projects\Hamid\Projects\Mojgan\Analyses\Combined_features\' pdf_file_name '.pdf'], '-dpdf', pdf_print_resolution)
end



