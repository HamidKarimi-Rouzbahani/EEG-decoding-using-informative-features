 function [signal] = DatasetLoading_cambridge(dataset,subject)
% dataset=1;
% subject=1;
if dataset==1
    %% Dataset 1, mine
    raster_file_directory_name = '/group/woolgar-lab/projects/Hamid/Projects/Mojgan/EEG_data/';
    if subject==1
        save_prefix_name = 'MahdiBorjian';
        session_date='20151124';
    elseif subject==2
        save_prefix_name = 'HamidGolmohammadi';
        session_date='20151201';
    elseif subject==3
        save_prefix_name = 'FarzadAghayeezade';
        session_date='20151201';
    elseif subject==4
        save_prefix_name = 'PouyaAhmadvand';
        session_date='20151208';
    elseif subject==5
        save_prefix_name = 'AliZarandi';
        session_date='20151215';
    elseif subject==6
        save_prefix_name = 'Khanum1';
        session_date='20151216';
    elseif subject==7
        save_prefix_name = 'Sub12';
        session_date='20160112';
    elseif subject==8
        save_prefix_name = 'Sub13';
        session_date='20160112';
    elseif subject==9
        save_prefix_name = 'Sub17';
        session_date='20160113';
    elseif subject==10
        save_prefix_name = 'Sub19';
        session_date='20160117';
    end
    
    raster_file_dir = dir([raster_file_directory_name,save_prefix_name,session_date,'*.mat']);
    channels=[1:length(raster_file_dir)];
    data=[];
    c=0;
    for chan=[1 12 23 26:31 2:11 13:22 24 25]
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
    raster_file_directory_name = '/group/woolgar-lab/projects/Hamid/Projects/Mojgan/Vahab dataset/';
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
    data_labels=load(['/group/woolgar-lab/projects/Hamid/Projects/Mojgan/StanfordDataset/S',num2str(subject),'_a1_2/CH',num2str(ch),'.mat'],'raster_labels');
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
        X=load(['/group/woolgar-lab/projects/Hamid/Projects/Mojgan/StanfordDataset/S',num2str(subject),'_a1_2/CH',num2str(ch),'.mat'],'raster_data');
        for tr=1:size(X.raster_data,1)
            X.raster_data(tr,:)=X.raster_data(tr,:)-nanmean(X.raster_data(tr,301:500));
        end
        for cat=1:6
            signal(ch,:,cat,1:length(trials_inds(cat,:)))=X.raster_data(trials_inds(cat,:),301:end)';
        end
    end
end
end
