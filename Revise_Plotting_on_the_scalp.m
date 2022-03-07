clc;
clear all;
close all;
Dataset=3; %Mine Vhab Stfd
Ranks_data=0; % ranks or values

if Dataset>2
    channels=128;
else
    channels=31;
end
Datasets={'Mine','Vhab','Stfd'};

for Subject=1:10
    load(['Revise_corrected_Dec_DS_',Datasets{Dataset},'_Band_Theta_Wind_sliding_Subject_',num2str(Subject),'_Wave_Camb.mat'],'accuracy','deltas');
    
    delta=squeeze(nanmean(deltas(2,:,:,:),2));
    [sorted_values,sorted_ranks]=sort(delta');
    for windoww=1:231
        times(windoww)=((windoww-1)*5+1+(windoww-1)*5+50)./2;
        for ch=1:channels
            inds=find(sorted_ranks(:,windoww)==ch);
            channel_ranks(ch,windoww,Subject)=inds;
        end
    end
    channel_values(:,:,Subject)=delta';
end
if Ranks_data==1
    X=squeeze(nanmean(channel_ranks,3))';
else
    X=squeeze(nanmean(channel_values,3))';
end
save('X_data.mat','X')
eeglab
