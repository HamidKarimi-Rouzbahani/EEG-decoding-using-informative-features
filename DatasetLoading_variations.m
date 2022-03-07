function [signal] = DatasetLoading_variations(dataset,subject,variation)

if dataset==1
    %% Dataset 1, mine
    raster_file_directory_name = 'D:\Hamid\Mojgan_analyses\My dataset\';
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

    
    if variation==1
        var_char='l';
    elseif variation==2
        var_char='p';        
    elseif variation==3
        var_char='s';        
    elseif variation==4
        var_char='t';       
    end
    
    c=0;
    for chan=[1 12 23 26:31 2:11 13:22 24 25]
        load([raster_file_directory_name raster_file_dir(chan).name]);
        c=c+1;
        data(:,:,c)=raster_data;
    end
    d=0;
    for j=1:size(raster_labels.combined_ID_subVariation,2)
        if strcmp(raster_labels.combined_ID_subVariation{1,j}(1,end-2),var_char)
            d=d+1;
            indxes(d)=j;
            chosen_var{1,d}=raster_labels.combined_ID_subVariation{1,j}(1,end-2:end);
        end
    end
    labels=chosen_var;
    data=data(indxes,:,:);
    [Categories,~,indices]=unique(labels);
    activity=cell(31,length(Categories));
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
    signal = nan.*ones(31,1200,length(Categories),fix(144./length(Categories)));    
    for ch=1:31
        for cat=1:length(Categories)
            for trial=1:size(activity{ch,cat},1)
                signal(ch,:,cat,trial)=activity{ch,cat}(trial,:);
            end
        end
    end
end
