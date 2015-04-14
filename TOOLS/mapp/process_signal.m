function process_signal(probe)
    global strct_config
    global strct_probes

    %%%Originalsignal beschneiden und reduzieren
    start=strct_config.firstSample;
    ende=strct_config.lastSample;
    strct_probes(probe).processedDataX=strct_probes(probe).rawDataX(start:ende);
    strct_probes(probe).processedDataY=strct_probes(probe).rawDataY(start:ende);
    strct_probes(probe).samples=length(strct_probes(probe).processedDataY);

    if strct_config.L_ref_flag==1

        if probe==length(strct_config.markiert)
            fprintf('Korrektur der Länge und somit der dimensionsbehafteten Frequenz von L=%g auf L=%g\n',strct_config.L_ref_config,strct_config.global_L_ref);
        end
        %Sr=f*L/U : Wenn die Daten vergleichbar sein sollen muss die zu
        %Dimensionierung benutzte L Größe zunächst herausgerechnet und
        %auf das globale L_ref korrigiert werden.
        L_ref_corrector=strct_config.L_ref_config/strct_config.global_L_ref;
        strct_probes(probe).processedDataX=strct_probes(probe).processedDataX/L_ref_corrector;

    end
    
  
    %Y-Originaldaten werden um den Mittelwert korrigiert, damit diese um 0
    %oszillieren
    if strct_config.avrg==1 || strct_config.co_frequency ~=0
        strct_probes(probe).processedDataY=strct_probes(probe).processedDataY-mean(strct_probes(probe).processedDataY);
    end
    
    %Sofern ausgewählt wird ein Highpassfilter auf die Daten angewendet
    if strct_config.co_frequency ~=0
        Fs=1/(strct_probes(probe).processedDataX(2)-strct_probes(probe).processedDataX(1));
        Fc=strct_config.co_frequency/(Fs/2); %Cutoff frequency
        [b,a]= butter(4,Fc,'high');
        tempData=filtfilt(b,a,strct_probes(probe).processedDataY);
        strct_probes(probe).processedDataY =tempData;
    end
    
    %Sofern ausgewählt wird das hanning-Fenster auf die Daten angewendet
    if strct_config.hanning==1
       % Using windows function
       w = hanning(strct_probes(probe).samples) ;
       CG=2; %Correction factor of Hanning window (see Paper Schmid2009)
       strct_probes(probe).processedDataY = w.*strct_probes(probe).processedDataY*CG ;
    end
    
end
