clc;
clear all;
% close all;
Dataset=1;

Subjects=[5];
bands=[2]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1];

Fs=1000;

%% These features provide 1 number for each trial and each channel
for band=bands   % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
    
    if band==2
        lowband=0.5;
        highband=4;
    elseif band==3
        lowband=4;
        highband=8;
    elseif band==4
        lowband=8;
        highband=12;
    elseif band==5
        lowband=12;
        highband=16;
    elseif band==6
        lowband=25;
        highband=200;
    end
    
    
    for Subject=Subjects
        signal = DatasetLoading(Dataset,Subject); % channel, time, cat, trial
        
        for windoww=Windows
            
            windows{1}=[201:1200]; % Post-stim
            windows{2}=[1:200];    % Pre-stim
            windows{3}=[200:300];  
            windows{4}=[301:400];
            windows{5}=[401:500];
            windows{6}=[501:600];
            windows{7}=[601:700];
            windows{8}=[701:800];
            
            wind=windows{windoww};
            
            
            %   for feature=[1:9 11:30 32:35] % features when evaluating
            %   broad-band whole-time windows: set band=1 and windoww=1
            
            %   for feature=[1:9 11:20 27:30 32:35] % features when evaluating freq bands
            
            %   for feature=[2:9 11:13 18:30 32 34:35] % features when evaluating time windows
            
            accuracy = nan*ones(35,nchoosek(size(signal,3),2),10);
            for feature=[34]
%               for feature=[2:9 11:13 18:30 32 34:35] % features when evaluating time windows
                
                if feature==33
                    net= load ('imagenet-caffe-alex.mat');
                end
                for channel = 1:size(signal,1)
                    for category = 1:size(signal,3)
                        for trial = 1:size(signal,4)
                            
                            if feature==1
                                trial_data=signal(channel,1:200,category,trial);
                            else
                                trial_data=signal(channel,wind,category,trial);
                            end
                            if band>1
                                trial_data=eegfilt(trial_data,Fs,lowband,highband,0,floor(length(trial_data)/3),0,'fir1');
                            end
                            
                            if feature<28 && channel==1 && category==1 && trial==1 % single-valued
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4));
                            elseif feature==28 && channel==1 && category==1 && trial==1 % Wavelet
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),length(wind)+20);
                            elseif (feature==29 || feature==30) && channel==1 && category==1 && trial==1 % Hilbert
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),length(wind));
                            elseif feature==31 && channel==1 && category==1 && trial==1 % PAC
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),size(signal,1)-1);
                            elseif feature==32 && channel==1 && category==1 && trial==1 % ICC
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),size(signal,1));
                            elseif feature==33 && channel==1 && category==1 && trial==1 % CNN
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),1000);
                            elseif feature==34 && channel==1 && category==1 && trial==1 % Samples
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),length(wind));
                            elseif feature==35 && Dataset<3 && channel==1 && category==1 && trial==1 % Autocorr
                                feature_chosen=nan.*ones(size(signal,1),size(signal,3),size(signal,4),size(signal,1));
                            end
                            
                            if sum(isnan(trial_data))>0
                                feature_chosen(channel,category,trial)=nan;
                            else
                                if feature==1
                                    %%            Baseline
                                    feature_chosen(channel,category,trial) = mean(trial_data);
                                elseif feature==2
                                    %%            Time features: Sigal mean
                                    feature_chosen(channel,category,trial) = mean(trial_data);
                                    
                                elseif feature==3
                                    %% Signal median
                                    feature_chosen(channel,category,trial) = median(trial_data);
                                    
                                elseif feature==4
                                    %% Signal variance
                                    feature_chosen(channel,category,trial) = var(trial_data);
                                    
                                elseif feature==5
                                    %% Signal skewness
                                    feature_chosen(channel,category,trial) = skewness(trial_data);
                                    
                                elseif feature==6
                                    %% Signal Kurtosis
                                    feature_chosen(channel,category,trial) = kurtosis(trial_data);
                                    
                                elseif feature==7
                                    %%                 LZ complexity
                                    threshold = median(trial_data);
                                    %         threshold = mean(data);
                                    trial_data(trial_data>= threshold)=1;
                                    trial_data(trial_data< threshold)=0;
                                    [feature_chosen(channel,category,trial),~] = calc_lz_complexity(trial_data, 'exhaustive', 1);
                                    
                                elseif feature==8
                                    %%                  Higuchi fractal dimension
                                    maxtime = length(trial_data);
                                    %                             Kmax = floor(maxtime./2);
                                    Kmax = 10;
                                    feature_chosen(channel,category,trial) = Higuchi_FD(trial_data,Kmax);
                                    %                 % second implementation
                                    %                 [HFD2(channel,category,trial),~,~,~]  = HFD_LCALC(data);
                                    %                 % thrid implementation
                                    %                 HFD3(channel,category,trial) = hfd(data,Kmax);
                                    
                                elseif feature==9
                                    %%                 Katz fractal dimensions
                                    feature_chosen(channel,category,trial) = Katz_FD(trial_data);
                                    
                                elseif feature==10
                                    %%                  Lyapunov exponent (largest LLE)
                                    [feature_chosen(channel,category,trial),~] = lyaprosen(trial_data,0,0);
                                    
                                elseif feature==11
                                    %%                  Hurst Exponent
                                    feature_chosen(channel,category,trial) = estimate_hurst_exponent(trial_data);
                                    %                 HE2(channel,category,trial) = genhurst(data);
                                    %                 HE3(channel,category,trial) = hurstCC(data);
                                    
                                elseif feature==12
                                    %%                  Sample entropy
                                    feature_chosen(channel,category,trial) = entropy (trial_data);
                                    %                 Ent2(channel,category,trial) = SampEn (2,0.2.* std(data),data,1);
                                    
                                elseif feature==13
                                    %%                  Approximate Entropy
                                    feature_chosen(channel,category,trial) = ApEn (2,0.2.* std(trial_data),trial_data,1);
                                    %                 Ent4(channel,category,trial) = approx_entropy(2,0.2.* std(data),data);
                                    
                                elseif feature==14
                                    %%                  ERP components P1/C1/P100, N1/N170, P200/P2a and P2b, chosen by visual inspection of PO3 channel, windows of components are 40, 60, 70 and 80 ms respectively.
                                    % each trial can have 1 to 4 components
                                    % P1
                                    feature_chosen(channel,category,trial)=nanmean(trial_data([80:120]));
                                    
                                elseif feature==15
                                    %% N1
                                    feature_chosen(channel,category,trial)=nanmean(trial_data([120:200]));
                                    
                                elseif feature==16
                                    %% P2a
                                    feature_chosen(channel,category,trial)=nanmean(trial_data([150:220]));
                                    
                                elseif feature==17
                                    %% P2b
                                    feature_chosen(channel,category,trial)=nanmean(trial_data([200:275]));
                                    
                                elseif feature==18
                                    %%              Within-trial correlation
                                    numLags=size(signal,1)+1;
                                    [acf,lags,~] =autocorr(trial_data,numLags);
                                    feature_chosen(channel,category,trial)= mean(acf(2:end));
                                    
                                elseif feature==19
                                    %%              Hjorth complexity
                                    % this finds spread of the spectrum and represents the change in frequency
                                    % Hcomplexity
                                    step_size=1./Fs;
                                    data_prime=(diff(trial_data)./step_size);
                                    data_second=(diff(data_prime)./step_size);
                                    feature_chosen(channel,category,trial)=(std(data_second).*std(trial_data))./(std(trial_data)).^2;
                                    
                                elseif feature==20
                                    %% Hmobility
                                    step_size=1./Fs;
                                    data_prime=(diff(trial_data)./step_size);
                                    data_second=(diff(data_prime)./step_size);
                                    feature_chosen(channel,category,trial)=std(data_prime)./std(trial_data);
                                    
                                elseif feature==21
                                    %% Mean Freq
                                    feature_chosen(channel,category,trial) = meanfreq(trial_data,Fs);
                                    
                                elseif feature==22
                                    %% Median Freq
                                    feature_chosen(channel,category,trial) = medfreq(trial_data,Fs);
                                    
                                elseif feature==23
                                    %% Average Freq
                                    zeroscount=0;
                                    for i=2:length(trial_data)
                                        if (trial_data(i)>0 && trial_data(i-1)<0) || (trial_data(i)<0 && trial_data(i-1)>0)
                                            zeroscount=zeroscount+1;
                                        end
                                    end
                                    feature_chosen(channel,category,trial) =zeroscount.*(length(trial_data)./Fs);
                                    
                                elseif feature==24
                                    %% Spectral edge frequency 95%
                                    if var(trial_data)==0
                                        feature_chosen(channel,category,trial) =0;
                                    else
                                        Fourier = fft(trial_data)/length(trial_data);
                                        Fouriers = (abs(Fourier));                           % Spectrum
                                        Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;        % Frequency Vector
                                        Iv = 1:length(Fv);                                      % Index Vector
                                        IntSpectrum = cumtrapz(Fv, Fouriers(Iv));               % Numeric Integration
                                        feature_chosen(channel,category,trial) = interp1(IntSpectrum, Fv, 0.95*IntSpectrum(end), 'linear');    % Interploate To Find ‘SEF’
                                    end
                                elseif feature==25
                                    %% Power at Median Freq
                                    amp = 2*abs(fft(trial_data))/length(trial_data);
                                    phs = angle(fft(trial_data));
                                    Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;  % Frequency Vector
                                    Iv = 1:length(Fv);                                      % Index Vector
                                    fr_des = medfreq(trial_data,Fs);                        % Desired Frequency
                                    ampv = amp(Iv);                                         % Trim To Length Of ‘Fv’
                                    phsv = phs(Iv);                                         % Trim To Length Of ‘Fv’
                                    ap = [ampv(:) phsv(:)];                                 % Amplitude & Phase Matrix
                                    ap_des = interp1(Fv(:), ap, fr_des, 'linear');
                                    feature_chosen(channel,category,trial) =ap_des(1);
                                    
                                elseif feature==26
                                    %% Phase at Median Freq
                                    amp = 2*abs(fft(trial_data))/length(trial_data);
                                    phs = angle(fft(trial_data));
                                    Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;  % Frequency Vector
                                    Iv = 1:length(Fv);                                      % Index Vector
                                    fr_des = medfreq(trial_data,Fs);                        % Desired Frequency
                                    ampv = amp(Iv);                                         % Trim To Length Of ‘Fv’
                                    phsv = phs(Iv);                                         % Trim To Length Of ‘Fv’
                                    ap = [ampv(:) phsv(:)];                                 % Amplitude & Phase Matrix
                                    ap_des = interp1(Fv(:), ap, fr_des, 'linear');
                                    feature_chosen(channel,category,trial) =ap_des(2);
                                    
                                elseif feature==27
                                    %% Signal power
                                    feature_chosen(channel,category,trial)=bandpower(trial_data);
                                    
                                elseif feature==28
                                    %% Wavelet transform:
                                    [c,l] = wavedec(trial_data,5,'sym2');
                                    [ca5] = appcoef(c,l,'sym2',5);
                                    [cd1,cd2,cd3,cd4,cd5] = detcoef(c,l,[1 2 3 4 5]);
                                    feature_chosen(channel,category,trial,1:length([ca5,cd1,cd2,cd3,cd4,cd5]))=[ca5,cd1,cd2,cd3,cd4,cd5];
                                    
                                elseif feature==29
                                    %% Hilbert transform amplitude
                                    Hilb = hilbert(trial_data);
                                    feature_chosen(channel,category,trial,:)=abs(Hilb);
                                    
                                elseif feature==30
                                    %% Hilbert transform phase
                                    Hilb = hilbert(trial_data);
                                    feature_chosen(channel,category,trial,:)=angle(Hilb);
                                    
                                elseif feature==31
                                    %% Phase-Amplitude Coupling
                                    data.trial{1,1}= trial_data;
                                    data.time{1,1}= [1:length(trial_data)]./Fs;
                                    data.trialinfo=[100]';
                                    data.label{1,1}='SampleData';
                                    toi=[0.001 1.0]; % time of interest
                                    phase=[0.5 12];   % phase(1):2.5:phase(2)
                                    ampl=[24 120];   % amp(1):19:amp(2)
                                    diag = 'no'; %'yes' or 'no' to turn on or off diagrams during computation
                                    surrogates = 'no'; %'yes' or 'no' to turn on or off surrogates during computation
                                    approach = 'tort';%,'ozkort','canolty','PLV';
                                    [MI_matrix_raw,~] = calc_MI(data,toi,phase,ampl,diag,surrogates,approach);
                                    feature_chosen(channel,category,trial,:)=reshape(MI_matrix_raw(1:6,1:5),[30 1]);
                                    
                                elseif feature==32
                                    %% Inter-channel correlation 31 feature per trial
                                    ICC=zeros(size(signal,1),1);
                                    for ch2=1:size(signal,1)
                                        trial_data2=signal(ch2,wind,category,trial);
                                        if band>1
                                            trial_data2=eegfilt(trial_data2,Fs,lowband,highband,0,floor(length(trial_data2)/3),0,'fir1');
                                        end
                                        ICC(ch2,1)=corr(trial_data',trial_data2');
                                    end
                                    feature_chosen(channel,category,trial,:)=ICC;
                                    
                                elseif feature==33
                                    %% CNN
                                    trial_data=(trial_data+abs(min(trial_data)));
                                    trial_data=(trial_data).*(255./max(trial_data));
                                    
                                    im_ = single(uint8(trial_data));
                                    im_ = repmat(im_,[1 ceil(227*227./length(im_))]);
                                    im_=im_(1,1:227*227);
                                    im_=reshape(im_,[227 227]);
                                    im_ = imresize(im_,net.meta.normalization.imageSize(1:2));
                                    im_ = im_ - net.meta.normalization.averageImage;
                                    
                                    res = vl_simplenn(net, im_);
                                    feature_chosen(channel,category,trial,:)=squeeze(res.x);
                                    
                                elseif feature==34
                                    %% Signal samples
                                    feature_chosen(channel,category,trial,:)=trial_data;
                                    
                                elseif feature==35
                                    %% Within-trial correlation
                                    numLags=size(signal,1);
                                    [acf,lags,~] =autocorr(trial_data,numLags);
                                    feature_chosen(channel,category,trial,:)= acf(2:end);
                                end
                            end
                        end
                    end
                end
                ccc
                %% Classification
                Xt=[];
                Y=[];
                for cat=1:size(signal,3)
                    isnotnan=~isnan(feature_chosen(:,cat,:,:));
                    tmp=squeeze(feature_chosen(:,cat,isnotnan(1,1,:,1),:));
                    tt=[];
                    for trial=1:sum(isnotnan(1,1,:,1))
                        tt=horzcat(tt,reshape(tmp(:,trial,:),[size(tmp,1)*size(tmp,3) 1]));
                    end
                    Xt=horzcat(Xt,tt);
                    Y=horzcat(Y,cat*ones(1,sum(isnotnan(1,1,:,1))));
                end
                X=Xt';
                Y=Y';

                if size(X,2)>size(signal,1)
                    coeff = pca(X');
                    X=coeff(:,1:size(signal,1));
                end
                
                clearvars -except Windows Subjects bands lowband highband band windows windoww wind Dataset Fs True_Predicted_labels accuracy Subject signal feature X Y  net
                
                folds=10;
                combinations=nchoosek(unique(Y),2);
                for combination=1:size(combinations,1)
                    c=0;
                    for classes=combinations(combination,:)
                        for counter=1:length(Y)
                            if Y(counter)==classes
                                c=c+1;
                                YY(c)=Y(counter);
                                XX(c,:)=X(counter,:);
                            end
                        end
                    end
                    Xready=XX;
                    Yready=YY;
                    
                    clearvars YY XX
                    ccc
                    Classifier_Model = fitcdiscr(Xready,Yready,'DiscrimType','pseudoLinear');
                    cvmodel = crossval(Classifier_Model);
                    L = kfoldLoss(cvmodel);
                    accuracy(feature,combination)=1-L;
                    
                    [band Subject windoww feature combination]
                end
                Datasets={'Mine','Vhab','Stfd'};
                Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
                save(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
            end
        end
    end
end

%% plotting
clc;
clear all;
close all;
bands=[1:6]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:8];
Datasets={'Mine','Vhab','Stfd'};
Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
band=1;
windoww=1;
Dataset=2;
for windoww=Windows
    for Subject=1:10
        load(['Corrected_Dec_DS_',Datasets{Dataset},'_Band_',Bands{band},'_Wind_',num2str(windoww),'_Subject_',num2str(Subject),'.mat'],'accuracy');
        accuracies(:,Subject)=nanmean(nanmean(accuracy,2),3);
    end

plot(nanmean(accuracies,2))
hold on;
end
