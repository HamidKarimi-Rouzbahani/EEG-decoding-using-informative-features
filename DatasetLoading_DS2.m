function [signal] = DatasetLoading_DS2(subject,cued)

if cued==1
    cts=1:4;
else
    cts=5:8;
end

ch=1;
data_labels=load(['D:\Hamid\Mojgan_analyses\Vahab_complete_data\Subj',num2str(subject),'\Ch (',num2str(ch),').mat'],'raster_labels');
load(['D:\Hamid\Mojgan_analyses\Vahab_complete_data\Subj',num2str(subject),'\Ch (',num2str(ch),').mat'],'raster_data');
labels=unique(data_labels.raster_labels.stimulus_CatCue);
cat=0;
for cats=cts
    cat=cat+1;
    c=0;
    for trial=1:size(raster_data,1)
        if strcmp(data_labels.raster_labels.stimulus_CatCue{1,trial},labels{1,cats})
            c=c+1;
            trials_inds(cat,c)=trial;
        end
    end
end
trials_inds(trials_inds==0)=nan;


signal = nan.*ones(31,1200,4,144);
for ch=1:31
    X=load(['D:\Hamid\Mojgan_analyses\Vahab_complete_data\Subj',num2str(subject),'\Ch (',num2str(ch),').mat'],'raster_data');
    for tr=1:size(X.raster_data,1)
        X.raster_data(tr,:)=X.raster_data(tr,:)-nanmean(X.raster_data(tr,1301:1500));
    end
    cat=0;
    for cats=cts
        cat=cat+1;
        signal(ch,:,cat,1:length(trials_inds(cat,~isnan(trials_inds(cat,:)))))=X.raster_data(trials_inds(cat,~isnan(trials_inds(cat,:))),1301:2500)';
    end
end
end
