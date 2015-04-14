%Maximum Entropy Method
%14.5.2014
% Berechnet die AR-Koeffizienten fï¿½r eine vorgegebene Model-Ordnung.
% Anschlieï¿½end wird der FPE nach Akaike berechnet und aus den Empfehlungen
% von [Ulrych] die beste Ordnung ausgesucht.
% Angezeigt wird dann die Leistungsdichte berechnet mit der MEM und FFT
%


%clc
%clear
%format longEng
%load('Testsignal','y')

function order_optimum = ar_order(y, model_order,probe)

%% Variablen
N = length(y);                                                                  % Signallï¿½nge
%model_order = 300;                                                              % Modelordnung kann hier eingestelle werden. !!!WICHTIG 2<M<N!!!
if model_order >= N || model_order < 2
    logbar(sprintf('Fehler: Ordnung muss zwischen 1 und %d liegen',N));
    order_optimum = 0;
else
logbar(sprintf('Berechne optimale AR - Ordnung für Probe %d zwischen 2 und %d ...',probe, model_order));
P = zeros(1,model_order);               
b = zeros((model_order+1),(N-model_order+1));       
a = zeros((model_order+1),(N-model_order+1));
alpha = zeros(model_order+1,N-model_order+1);                                   % AR Koeffizient - [M x k]
alpha(1:model_order+1) = -1;
%%Initialisierung von M=1
for j=1:1:(N-1)    
    A(1,j) = y(1,j) * y(1,j+1);    
    B(1,j) = (y(1,j))^2 + (y(1,j+1))^2;    
end
ZA = 2*sum(A);
NE = sum(B);
alpha(2,2)=ZA/NE;
clear ZA NE A B
A = 1/(2*(N-1));
for j=1:1:(N-1)    
    B(1,j) = (y(1,j+1)-alpha(2,2)*y(1,j))^2;    
    C(1,j) = (y(1,j)- alpha(2,2)*y(1,j+1))^2;    
end
SUMMEB = sum(B);
SUMMEC = sum(C);
P(1,2) = A*(SUMMEB+SUMMEC);
clear A B C SUMMEB SUMMEC
R(1,1) = P(1,2)/(1-alpha(2,2)^2);                                           %R_0
R(1,2) = R(1,1)*alpha(2,2);                                                 %R_1

%% Schrittweise Erhï¿½hung von M, berechnung der AR-Koeffizienten und P(M+1)
for M=2:1:model_order
    %% Berechnung von alpha(M,M)
    %logbar(sprintf('Berechne FPE für Ordnung %d',M));
    for j=1:1:(N-M)
        for k=1:1:(M+1)            
            B(1,k) = alpha(M,k)*y(1,j+k-1);
        end    
        b(M+1,j+1) = sum(B);
        clear B
        for k=0:1:M
            A(1,k+1) = alpha(M,M-k+1)*y(1,j+k);
        end
        a(M+1,j+1) = sum(A);  
        clear A
    end
    for j=1:1:(N-M)
        ze(1,j) = b(M+1,j+1)*a(M+1,j+1);        
        ne(1,j) = (b(M+1,j+1))^2+(a(M+1,j+1))^2;
    end    
    ZE = 2*sum(ze);
    NE = sum(ne);    
    alpha(M+1,M+1) = ZE/NE;                                                 %alpha(M,M)
    clear ze ne    
    %% Berechnung der restlichen alpha-Werte [=> alpha(M,k)]    
    for k=1:1:(M-1)
        alpha(M+1,k+1) = alpha(M,k+1) - alpha(M+1,M+1)*alpha(M,M+1-k); 
    end
    
    
    for j=1:1:(N-M)
    B(1,j) = (b(M+1,j+1)-alpha(M+1,M+1)*a(M+1,j+1))^2 + (a(M+1,j+1) - alpha(M+1,M+1)*b(M+1,j+1))^2;
    end
    C = sum(B);
    A = 1/(2*(N-M));
    PM_1(1,M-1) = A*C;
    clear A B C
    
    FPE_LINKS(1,M-1) = (N+M+1)/(N-M-1);
    
    
    FPE(1,M-1)=PM_1(1,M-1)*FPE_LINKS(1,M-1);
    
    
    
end
[~,order_optimum] = min(FPE);
logbar(sprintf('FPE minimal bei Ordnung %d',order_optimum));
end

% figure
% subplot(2,3,1)
%     plot(y)
%     grid on
% subplot(2,3,2)
%     hold on
%     plot(2:M,FPE/max(FPE),'color','red')
%     plot(2:M,FPE_LINKS/max(FPE_LINKS),'color','green')
%     plot(2:M,PM_1/max(PM_1),'color','blue')
%     legend('FPE','Linker Therm','P_m+1')
%     xlim([2 model_order])
% subplot(2,3,3)
%     h_burg = spectrum.burg(order_optimum);    
%     hpsd_burg = psd(h_burg,y)%,'Fs',Fs);
%     plot(hpsd_burg)   
%     title('MEM')
% subplot(2,3,4)
%     L = N;
%     NFFT = 2^nextpow2(L);
%     Y = fft(y,NFFT)/L;
%     plot(2*abs(Y(1:NFFT/2+1)))
%     title('FFT')
%     grid on
% subplot(2,3,5)
%     nfft=length(y);
%     if mod(nfft,2)==0
%         size=nfft/2+1;
%     else
%         size=(nfft+1)/2;
%     end
%     fs=85833;
%         
%     pburgY=zeros(size);
%     pburgF = fs*linspace(0,1,size);
%     [pburgY,pburgF]=pburg(y,order_optimum,pburgF,fs);
%     plot(pburgF,pburgY)
%     axis([800 5e3 0 10])
%     axis 'auto y'
%     grid on
end







