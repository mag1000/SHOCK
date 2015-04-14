function AuswertungStarten(path)
    global CoordinateX
    global CoordinateY
    global importedData
    global samples
    global probes
    global data
    global flag_psdTKE
    global psd_Y
    global psd_f    
    
    global flag_FFT
    
    global flag_spectrogram
    global spectrogram_S
    global spectrogram_F
    global spectrogram_T
    global flag_periodogram

    global flag_pwelch
    
    global flag_pburg    
        
    global flag_13oktave
    
    
    
    importedData=importData(path);
    importedData(1,:)=importedData(1,:)-importedData(1,1);
    time=importedData(1,:);
    %delta_t=0.000045;
    %time=linspace(0,samples*delta_t,samples);
    Fs=1/((time(length(time))-time(1))/(samples-1));
    fprintf('Fs: %g\n',Fs);
    fprintf('Es werden %d samples eingelesen.\n',samples);
    
    
    %f1=7.5*sawtooth(time*7500*pi);
%     f1=13.3*sin(2*pi*2000*time);
%     f2=5.2*sin(2*pi*1000*time);
%     data=f1+f2+2*randn(size(time)); % Sinusoids plus noise
%     for index=1:probes
%     importedData(index+1,:)=data;
%     end
    
    if flag_psdTKE==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%   psdTKE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nfft=2^nextpow2(samples);
    %nfft=samples;
    data=importedData(2,:);
    [pxx_temp,f_temp] = periodogram(data,rectwin(length(data)),nfft,Fs);
    psd_Y=zeros(probes+1,length(pxx_temp));
    for index=1:probes
        % Using windows function
        window = hanning(samples,'periodic');
        %window = rectwin(samples);
        data=importedData(index+1,:);
        data=data-mean(data);
        CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
        %data = window'.*data*CG; 
        [psd_Y(index,:),psd_f] = periodogram(data,window,nfft,Fs,'power');
        if index<probes
            dA=((CoordinateX(index)-CoordinateX(index+1))^2+(CoordinateY(index)-CoordinateY(index+1))^2)^0.5;
        else
            dA=((CoordinateX(index)-CoordinateX(index-1))^2+(CoordinateY(index)-CoordinateY(index-1))^2)^0.5;
        end
        
        psd_Y(probes+1,:)=psd_Y(probes+1,:)+psd_Y(index,:)*dA;
    end
    dA=((CoordinateX(probes)-CoordinateX(1))^2+(CoordinateY(probes)-CoordinateY(1))^2)^0.5;
    fprintf('Die Gesamthöhe der Auswertung beträgt :%f\n',dA);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
if flag_periodogram==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%   periodogram   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nfft=2^nextpow2(samples/1);
    %nfft=samples;
    data=importedData(2,:);
    [pxx_temp,f_temp] = periodogram(data,rectwin(length(data)),nfft,Fs);
    psd_Y=zeros(probes+1,length(pxx_temp));
    for index=1:probes
        % Using windows function
        window = hamming(samples,'periodic');
        %window = rectwin(samples);
        data=importedData(index+1,:);
        data=data-mean(data);
        CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
        %data = window'.*data*CG; 
        [psd_Y(index,:),psd_f] = periodogram(data,window,nfft,Fs);
        if index<probes
            dA=((CoordinateX(index)-CoordinateX(index+1))^2+(CoordinateY(index)-CoordinateY(index+1))^2)^0.5;
        else
            dA=((CoordinateX(index)-CoordinateX(index-1))^2+(CoordinateY(index)-CoordinateY(index-1))^2)^0.5;
        end
        
        psd_Y(probes+1,:)=psd_Y(probes+1,:)+psd_Y(index,:)*dA;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end    


    if flag_FFT==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%   FFT   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fft geteilt durch Anzahl sample liefert die amplitudenverteilung mit der
    %Einheit des Eingangsignals
    
    nfft=pow2(nextpow2(samples));
    %nfft=samples;
    psd_Y=zeros(probes+1,nfft);
    
    for index=1:probes
        % Using windows function
        %w = hamming(samples);
        w = hamming(samples,'periodic');
        data=importedData(index+1,:);
        data=data-mean(data);
        CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
        data2 = w'.*data*CG; 
        %data2 = data; 
        psd_Y(index,:) = 2*abs(fft(data2,nfft))/samples;
        %psd aus fft
        %psd_Y(index,:)=10*log10(psd_Y(index,:).*psd_Y(index,:)./length(psd_Y(index,:)));
        
        if index<probes
            dA=((CoordinateX(index)-CoordinateX(index+1))^2+(CoordinateY(index)-CoordinateY(index+1))^2)^0.5;
        else
            dA=((CoordinateX(index)-CoordinateX(index-1))^2+(CoordinateY(index)-CoordinateY(index-1))^2)^0.5;
        end
        
        psd_Y(probes+1,:)=psd_Y(probes+1,:)+psd_Y(index,:);%*dA;
    end
    
    psd_f = 1/(time(2)-time(1))*linspace(0,1,nfft);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    if flag_spectrogram==1
    %% %%%%%%%%%%%%%%%%%%%%%%%   Spectogram   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    US=2000;    %usedSamples
    CF=round(samples/US)-1;     %coarseFactor
    windowsize=256;
    window=hanning(windowsize);
    nfft=256;
    noverlap=windowsize-1;
    time2=time(1:CF:end);
    fs=1/(time2(2)-time2(1));
    
    for i=1:US
        data(i)=importedData(2,i*CF);
    end
    data=data-mean(data);
    [S,F,T] =spectrogram(data,window,noverlap,nfft,fs);
        
    spectrogram_S=zeros(length(F),length(T),probes);
    spectrogram_T=zeros(probes,length(T));
    spectrogram_F=zeros(length(F),probes);
    
    for index=1:probes
        for i=1:US
            data(i)=importedData(index+1,i*CF);
        end
        data=data-mean(data);
        [spectrogram_S(:,:,index),spectrogram_F(:,index),spectrogram_T(index,:)] =...
            spectrogram(data,window,noverlap,nfft,fs);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    end
    
    
    if flag_pwelch==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%   PWelch   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    fs=1/(time(2)-time(1));
    
    windowsize=pow2(nextpow2(samples/3));
    overlap=windowsize/2;

    fprintf('PWelch: Windowsize: %d, OverlappingPoints:%d\n',windowsize,overlap);
    window=hamming(windowsize);
    CG=2.;
    data=importedData(2,:);
    [pxx_temp,pxx_f] = (pwelch(data,window,overlap,[],fs,'psd'));
    psd_Y=zeros(probes+1,length(pxx_temp));
    for index=1:probes
        data=CG*importedData(index+1,:);
        data=data-mean(data);
        [psd_Y(index,:),psd_f] = pwelch(data,window,overlap,[],fs,'psd');
        
        psd_Y(index,:)=(psd_Y(index,:));
        %psd_Y(index,:)=db2mag(psd_Y(index,:));
        
        if index<probes
            dA=((CoordinateX(index)-CoordinateX(index+1))^2+(CoordinateY(index)-CoordinateY(index+1))^2)^0.5;
        else
            dA=((CoordinateX(index)-CoordinateX(index-1))^2+(CoordinateY(index)-CoordinateY(index-1))^2)^0.5;
        end
        
        psd_Y(probes+1,:)=psd_Y(probes+1,:)+psd_Y(index,:)/probes;
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end    
    
    if flag_13oktave==1
    %% %%%%%%%%%%%%%%%%%%%%%%%   1/3 Oktave SPL   %%%%%%%%%%%%%%%%%%%%%%%%%%
    BandsPerOctave = 1;
    N = 6;           % Filter Order
    F0 = 1000;       % Center Frequency (Hz)
    filter1 = fdesign.octave(BandsPerOctave,'Class 1','N,F0',N,F0,Fs);
    F0 = validfrequencies(filter1);
    Nfc = length(F0);
    for i=1:Nfc,
        filter1.F0 = F0(i);
        Hd(i) = design(filter1,'butter');
    end

    filter1.BandsPerOctave = 3;
    filter1.FilterOrder = 8;
    F0 = validfrequencies(filter1);
    Nfc = length(F0);
    for i=1:Nfc,
        filter1.F0 = F0(i);
        Hd3(i) = design(filter1,'butter');
    end    
    
    [Syyw, fPyyw] = pwelch(data,hamming(64),[],[],Fs,'power');
    yw = zeros(NFFT,Nfc);
    Pyyw = zeros(1,Nfc);
    for i=1:Nfc,
        yw(:,i) = filter(Hd3(i),data);
        [pPwr, pFreq] = pwelch(yw(:,i),hamming(64),[],[],Fs);
        Pyyw(i) = bandpower(pPwr, pFreq,'psd');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    if flag_pburg==1
    %% %%%%%%%%%%%%%%%%%%%%%%%  PBurg   %%%%%%%%%%%%%%%%%%
    nfft=samples;
    if mod(nfft,2)==0
        size=nfft/2+1;
    else
        size=(nfft+1)/2;
    end
    
    fs=1/(time(2)-time(1));
    model_order = 400;
    
    psd_Y=zeros(probes,size);
    psd_f = fs*linspace(0,1,size);
    
    for index=1:probes
        data=importedData(index+1,:);
        data=data-mean(data);
        [psd_Y(index,:),psd_f]=pburg(data,model_order,nfft,fs);

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
  

end

