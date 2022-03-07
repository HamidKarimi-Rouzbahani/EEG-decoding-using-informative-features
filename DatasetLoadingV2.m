% function [signal] = DatasetLoading(dataset,subject)

dataset=1;
subject=1;
if dataset==1
    %% Dataset 1, mine
    raster_file_directory_name = 'C:\Dataset1\';
    
    raster_file_dir = dir([raster_file_directory_name,'Subj',num2str(subject,'%02.f'),'*.mat']);
    channels=[1:length(raster_file_dir)];
    data=[];
    c=0;
    for chan=[1:31]
        load([raster_file_directory_name raster_file_dir(chan).name]);
        labels=raster_labels.stimulus_ID;
        c=c+1;
        data(:,:,c)=raster_data;
    end
    
    [Categories,~,indices]=unique(labels);
    activity=cell(31,4);
    spans=[1300:2499];
    chans=0;
    for n=channels
        chans=chans+1;
        category=0;
        for c=1:length(Categories)
            category=category+1;
            target_trials=find(indices==c);
            for trials=1:length(target_trials)
                activity{chans,category}(trials,:)=data(target_trials(trials),spans,n)-squeeze(mean(data(target_trials(trials),[1300:1500],n),2));
            end
        end
    end
    %% Dataset reformatting
    signal = nan.*ones(31,1200,4,144);
    for ch=1:31
        for cat=1:4
            for trial=1:size(activity{ch,cat},1)
                signal(ch,:,cat,trial)=activity{ch,cat}(trial,:);
            end
        end
    end
elseif dataset==2
    raster_file_directory_name = 'C:\Dataset2\';
    load([raster_file_directory_name ['S',num2str(subject),'.mat']]);
    for ch=1:31
        for cat=1:4
            for tr=1:size(FormattedData,4)
                FormattedData(ch,:,cat,tr)=FormattedData(ch,:,cat,tr)-nanmean(FormattedData(ch,1301:1500,cat,tr));
            end
        end
    end    
    signal=FormattedData(:,1300:2499,:,:);
    
elseif dataset==3
    ch=1;
    data_labels=load(['C:\Dataset3\S',num2str(subject),'_a1_2\CH',num2str(ch),'.mat'],'raster_labels');
    labels=unique(data_labels.raster_labels.Stimulus_Cat);
    for cat=1:6
        c=0;
        for trial=1:864
            if strcmp(data_labels.raster_labels.Stimulus_Cat{1,trial},labels{1,cat})
                c=c+1;
                trials_inds(cat,c)=trial;
            end
        end
    end
    
    signal = nan.*ones(128,1200,6,144);
    for ch=1:128
        X=load(['C:\Dataset3\S',num2str(subject),'_a1_2\CH',num2str(ch),'.mat'],'raster_data');
        for tr=1:size(X.raster_data,1)
            X.raster_data(tr,:)=X.raster_data(tr,:)-nanmean(X.raster_data(tr,301:500));
        end
        for cat=1:6
            signal(ch,:,cat,1:length(trials_inds(cat,:)))=X.raster_data(trials_inds(cat,:),301:end)';
        end
    end
end
end
