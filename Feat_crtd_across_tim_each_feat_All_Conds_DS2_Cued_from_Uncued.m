clc;
clear all;
% close all;
Subjects=[1:5];
bands=[1]; % 1=broad, 2=delta, 3=theta, 4=alpha, 5=beta, 6=gamma
Windows=[1:53];
Fs=1000;
Dataset=2;

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
        
        accuracy = nan*ones(35,length(Windows),16);
        for windoww=Windows
            
            windows=[1:1200];
            step_size=20;
            window_span=50;
            winds=windows((windoww-1)*step_size+1:(windoww-1)*step_size+window_span);
            
            features_chosen=nan*ones(31,4,288,21);
            for feature=[2:9 11:13 18:27] % features when evaluating time windows
                Xtotal=ones(1,31);
                Ytotal=0;
                for cued=0:1
                    signal = DatasetLoading_DS2(Subject,cued); % channel, time, cat, trial
                    
                    %                 for feature=[2] % features when evaluating time windows
                    if feature==33
                        net= load ('imagenet-caffe-alex.mat');
                    end
                    for channel = 1:size(signal,1)
                        for category = 1:size(signal,3)
                            for trial = 1:size(signal,4)
                                
                                if feature==1
                                    trial_data=signal(channel,1:200,category,trial);
                                else
                                    trial_data=signal(channel,winds,category,trial);
                                end
                                
                                if band>1 && sum(isnan(trial_data))==0
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
                                        %% Baseline
                                        feature_chosen(channel,category,trial) = mean(trial_data);
                                    elseif feature==2
                                        %% Time features: Sigal mean
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
                                        numLags=ceil(length(trial_data)./2);
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
                                            feature_chosen(channel,category,trial) = interp1(IntSpectrum, Fv, 0.95*IntSpectrum(end), 'linear');    % Interploate To Find �SEF�
                                        end
                                    elseif feature==25
                                        %% Power at Median Freq
                                        amp = 2*abs(fft(trial_data))/length(trial_data);
                                        phs = angle(fft(trial_data));
                                        Fv = linspace(0, 1, fix(length(trial_data)/2)+1)*Fs/2;  % Frequency Vector
                                        Iv = 1:length(Fv);                                      % Index Vector
                                        fr_des = medfreq(trial_data,Fs);                        % Desired Frequency
                                        ampv = amp(Iv);                                         % Trim To Length Of �Fv�
                                        phsv = phs(Iv);                                         % Trim To Length Of �Fv�
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
                                        ampv = amp(Iv);                                         % Trim To Length Of �Fv�
                                        phsv = phs(Iv);                                         % Trim To Length Of �Fv�
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
                                            trial_data2=signal(ch2,winds,category,trial);
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
                    %% Classification
                    
                    for ch=1:size(signal,1)
                        Xt=1;
                        Y=1;
                        for cat=1:size(signal,3)
                            isnotnan=~isnan(feature_chosen(:,cat,:,:));
                            tmp=squeeze(feature_chosen(ch,cat,isnotnan(1,1,:,1),:));
                            Xt=vertcat(Xt,reshape(tmp,[size(tmp,1)*size(tmp,2),1]));
                            Y=vertcat(Y,(cat+cued*4).*ones(size(tmp,1)*size(tmp,2),1));
                        end
                        X(1:length(Xt),ch)=Xt;
                    end
                    
                    clearvars -except Ytotal Xtotal Windows cued Subjects bands lowband highband band windows windoww winds Dataset Fs True_Predicted_labels accuracy Subject signal feature X Y  net
                    Xtotal=vertcat(Xtotal,X(2:end,:));
                    Ytotal=vertcat(Ytotal,Y(2:end,:));
                end
                Xtotal=Xtotal(2:end,:);
                Ytotal=Ytotal(2:end);
                combinations=[1 5;1 6;1 7;1 8;2 5;2 6;2 7;2 8;3 5;3 6;3 7;3 8;4 5;4 6;4 7;4 8];
                for combination=1:size(combinations,1)
                    c=0;
                    for classes=combinations(combination,:)
                        for counter=1:length(Ytotal)
                            if Ytotal(counter,1)==classes
                                c=c+1;
                                YY(c)=Ytotal(counter,1);
                                XX(c,:)=Xtotal(counter,:);
                            end
                        end
                    end
                    Xready=XX;
                    Yready=YY;
                    clearvars YY XX
                    Classifier_Model = fitcdiscr(Xready,Yready,'DiscrimType','pseudoLinear');
                    accuracy(feature,windoww,combination)=1-kfoldLoss(crossval(Classifier_Model));
                    [Subject windoww feature combination]
                end
                clearvars Xtotal Ytotal
            end
        end
        Datasets={'Mine','Vhab','Stfd'};
        Bands={'Broad','Delta','Theta','Alpha','Betta','Gamma'};
        save(['Cortd_Dec_DS_',Datasets{Dataset},'_Cued_vs_Uncued_Band_',Bands{band},'_Wind_slid_Subject_',num2str(Subject),'.mat'],'accuracy');
    end
end