
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Signalanalyse  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statistics(handles,probe)

global strct_probes
global strct_config

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strct_config.methode==1

    %fft geteilt durch Anzahl sample liefert die amplitudenverteilung mit der
    %Einheit des Eingangsignals
    
%     if isempty(handles) == 1
%         fprintf('FFT von Probe %s\n', strct_probes(probe).name);
%     else
%         logbar(handles.status_bar,sprintf('FFT von Probe %s', strct_probes(probe).name));
%     end
    
    nfft=pow2(nextpow2(strct_probes(probe).samples));
    
    strct_probes(probe).psdY=zeros(nfft/2+1,1) ;
    strct_probes(probe).psdX=zeros(nfft/2+1,1) ;
    
            
    full_fft_y = 2*abs(fft(strct_probes(probe).processedDataY,nfft))/nfft ;
    strct_probes(probe).psdY=full_fft_y(1:nfft/2+1);
    full_fft_x=1/strct_probes(probe).dx*linspace(0,1,nfft) ;
    strct_probes(probe).psdX=full_fft_x(1:nfft/2+1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pwelch%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strct_config.methode==3

    %fft geteilt durch Anzahl sample liefert die amplitudenverteilung mit der
    %Einheit des Eingangsignals
%     if isempty(handles) == 1
%         fprintf('PWelch von Probe %s\n', strct_probes(probe).name);
%     else
%         logbar(handles.status_bar,sprintf('PWelch von Probe %s', strct_probes(probe).name));
%     end
    
    strct_probes(probe).samples=length(strct_probes(probe).processedDataY);
    nfft=pow2(nextpow2(strct_probes(probe).samples)) ;
    
    strct_probes(probe).psdY=zeros(nfft,1) ;
    strct_probes(probe).psdX=zeros(nfft,1) ;
    
    %Da die pwelch Funktion das Signal in Segmente zerlegt, ist die
    %Fensterfunktion (sofern angewand) zunächst wieder zu entfernen bevor
    %die Daten verarbeitet werden können
    if strct_config.hanning==1
       % Using windows function
       w = hanning(length(strct_probes(probe).processedDataY)) ;
       CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
       strct_probes(probe).processedDataY = strct_probes(probe).processedDataY./w/CG;
    end
    
    windowsize=pow2(nextpow2(strct_probes(probe).samples/strct_config.pwelch_windowsize));
    overlap=windowsize/2;
    if strct_config.hanning==1
        window=hanning(windowsize);
    else
        window=rectwin(windowsize);
    end
    [strct_probes(probe).psdY,strct_probes(probe).psdX] = pwelch(strct_probes(probe).processedDataY,window,overlap,[],1/strct_probes(probe).dx,'psd');
    
    %Umwandlung PSD zu Amplitude:
    dx2=(strct_probes(probe).psdX(2)-strct_probes(probe).psdX(1));
    strct_probes(probe).psdY = sqrt(dx2*strct_probes(probe).psdY);
 
    %HanningFenster wird wieder hinzugefügt
    if strct_config.hanning==1
       % Using windows function
       w = hanning(length(strct_probes(probe).processedDataY)) ;
       CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
       strct_probes(probe).processedDataY = w.*strct_probes(probe).processedDataY*CG ;       
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%pBurg%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strct_config.methode==5
%     if isempty(handles) == 1
%         fprintf('PBurg von Probe %s\n', strct_probes(probe).name);
%     else
%         logbar(handles.status_bar,sprintf('PBurg von Probe %s', strct_probes(probe).name));
%     end
        strct_probes(probe).samples=length(strct_probes(probe).processedDataY);    

        nfft=pow2(nextpow2(strct_probes(probe).samples)) ;
        
        if mod(nfft,2)==0
            size=nfft/2+1;
        else
            size=(nfft+1)/2;
        end
    
        fs=1/(strct_probes(probe).processedDataX(2)-strct_probes(probe).processedDataX(1));
        
        strct_probes(probe).psdY = zeros(size);
        strct_probes(probe).psdX = fs*linspace(0,1,size);
        
        model_order = 100;
        [strct_probes(probe).psdY,strct_probes(probe).psdX]=pburg(strct_probes(probe).processedDataY,model_order,nfft,fs);
%         for i=1:length(strct_probes(1).psdY)
%         strct_probes.psdY=sqrt(strct_probes.psdY)*sqrt(2);
%         end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sigma%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strct_config.methode==6

%     if isempty(handles) == 1
%         fprintf('Standard deviation von Probe %s\n', strct_probes(probe).name);
%     else
%         logbar(handles.status_bar,sprintf('Standard deviation von Probe %s', strct_probes(probe).name));
%     end
      
    strct_probes(probe).psdY=std(strct_probes(probe).processedDataY);
    strct_probes(probe).psdX=strct_probes(probe).CoordinateX;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

        
        