function [Behaviour_per_category,RT_per_category] = DatasetLoading_DS2_behaviour(subject,cued)
if cued==1
    cts=1:4;
else
    cts=5:8;
end

ch=1;
data_labels=load(['Z:\projects\Hamid\\Projects\Mojgan\Vahab_complete_data\Subj',num2str(subject),'\Ch (',num2str(ch),').mat'],'raster_labels');
load(['Z:\projects\Hamid\\Projects\Mojgan\Vahab_complete_data\Subj',num2str(subject),'\Ch (',num2str(ch),').mat'],'raster_data');
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

trial_result = nan.*ones(4,72);
for cat=1:4
    for tr=1:length(trials_inds(cat,~isnan(trials_inds(cat,:))))
        if strcmp(data_labels.raster_labels.stimulus_Answer{1,trials_inds(cat,tr)}(end-7),'-')
            trial_result(cat,tr)=1;
        else
            trial_result(cat,tr)=0;
        end
    end
end

Behaviour_per_category=nansum(trial_result').\sum(~isnan(trial_result)');
%% RT
cts=1:4;
load(['Z:\projects\Hamid\\Projects\Mojgan\Vahab_complete_data\CogModPsychoRes\Subject',num2str(subject),'.mat']);
cat=0;
trials_inds_cor=[];
for cats=cts
    cat=cat+1;
    c=0;
    for trial=1:size(raster_data,1)
        if strcmp(data_labels.raster_labels.stimulus_CatCue{1,trial},labels{1,cats}) && strcmp(data_labels.raster_labels.stimulus_Answer{1,trial}(end-7),'-')
            c=c+1;
            trials_inds_cor(cat,c)=trial;
        end
    end
end
trials_inds_cor(trials_inds_cor==0)=nan;


trial_RT = nan.*ones(4,72);
for cat=1:4
    for tr=1:length(trials_inds_cor(cat,~isnan(trials_inds_cor(cat,:))))
        trial_RT(cat,tr)=CogModPsychoRes{trials_inds_cor(cat,tr),6};
    end
end
RT_per_category=nanmean(trial_RT,2);

end
