% berechnet die zeitliche Lösung der Orr-Sommerfeldgleichung für alpha=real
% und schätzt die räumliche Lösung nach gaster ab

function [f_real,alpha_i] = LST_NACA(profile)

global strct_config strct_probes

format long; 

spatial_flag=0; %räumliche Analyse (Wielandt)
blasius_flag=0; %Blasius Profil verwenden
flag_found_instability=0;

if strct_config.L_ref_flag==1
    L_ref_DNS=strct_config.global_L_ref; %m
else
    L_ref_DNS=strct_config.L_ref_config; %m
end

U_ref_DNS=strct_config.U_ref_config; %m/s

if blasius_flag~=1
    
    fid=fopen(profile,'r');
    fgetl(fid);
    strct_probes(1).name=strrep(sscanf(fgetl(fid),'%*s %s %*s %*s %*s %*s'), 't="','');
    strct_probes(1).name=strrep(strct_probes(1).name, '",','');
    DNS_Daten = fscanf(fid, '%g %g %g %g %g %g', [6 inf]);
    n_DNS=length(DNS_Daten);
    
    DNS_Y=DNS_Daten(1,:);
    DNS_U=DNS_Daten(2,:);
    DNS_delta1=DNS_Daten(3,1);
    DNS_delta=DNS_Daten(4,1);
    Re_delta1_DNS=DNS_Daten(5,1);    
    DNS_U0=DNS_Daten(6,1);
end

for i=1:1:n_DNS
    DNS_Y(i)=DNS_Y(i)/DNS_delta1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Orr-Sommerfeld-Solver %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_GS                = 50;
n                   = 300;

nu                  = 1.5*10^(-5); % Initialisierung der kinematische 
                                   % Viskosit�t von Luft bei 20�C in
                                   % [m^2/s]: nu

%U_inf               = DNS_U_infty*Ma_DNS*sqrt(1.4*287*300); 
                                   % Initialisierung der ungest�rten 
                                   % Anstr�mgeschwindigkeit in [m/s]: U_inf
                          
L_ref               = DNS_delta1 * L_ref_DNS;


delta_delta1        = DNS_delta/DNS_delta1;   % Berechnung des Verh�ltnis von 
                                              % Grenzschichtdicke zu Verdr�ngungdicke: 
                                              % delta_delta1
if blasius_flag==1, delta_delta1 = 5.0/1.7208; end ;       

delta_y             = delta_delta1/n_GS;      % Berechnung der Schrittweite
                                              % zwischen zwei Gitterpunkten
                                              % in y-Richtung mit 30% der 
                                              % Punkte innerhalb der 
                                              % Grenzschicht: delta_y

y_norm              = zeros(1,n+2); % Definition des Vektors der 
                                    % dimensionslosen y-Koordinate: y_norm
                                    
y_norm(1)           = 0; % Zuweisen des ersten Wertes der dimensionslosen     

                         % y-Koordinate: y_norm
                         

% Berechnen des Vektors der dimensionslosen y-Koordinate: y_norm %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:(n+2)                                                    %
                                                                 %
    y_norm(i)           = y_norm(i-1)+delta_y;                   %
                                                                 %
end                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if blasius_flag==1, Re_delta1_DNS=1000; end; %Re_delta1 für Blasius

Re_delta1               = Re_delta1_DNS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha_start             = strct_config.alpha_start; 
                           % Initialisierung des kleinsten zu berechnenden 
                           % Wertes der Wellenzahl alpha      
                           
alpha_end               = strct_config.alpha_end; 
                           % Initialisierung des gr��ten zu berechnenden 
                           % Wertes der Wellenzahl alpha
                           
anzahl_alpha_schritte   = strct_config.alpha_steps;
                           % Initialisierung der Anzahl der Schritte 
                           % für alpha

alpha_schritt           = (alpha_end-alpha_start)/(anzahl_alpha_schritte-1); 
                            % Initialisierung der Schrittweite zwischen 
                            % zwei Wellenzahlen alpha                           
                           
iterations_alpha        = round((alpha_end-alpha_start)/alpha_schritt+1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
number_amplified_eigenvalues = 0;   % Initialisierung des Z�hlers f�r die Vektoren die
                                    % ausschlie�lich instabile Punkte enthalten:
                                    % inst_alpha_final, inst_Re_delt1_final, 
                                    % inst_re_omega_final     
                                    
omega_alt                    = 0.2; % Initialisierung des Vergleichswertes für Omega,
                                    % wenn mehrere Eigenvektoren für ein alpha existieren
                                    
l                            = 1;   % Initialisierung des Z�hlers der Wellenzahlen                                    
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0                  = zeros(1,n); % Definition des Vektors der 
                                  % Koeffizienten a0 f�r die
                                  % Koeffizientenmatrix A0

a2                  = zeros(1,n); % Definition des Vektors der 
                                  % Koeffizienten a2 f�r die
                                  % Koeffizientenmatrix A2     
                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
% Berechnung der Ableitungsmatrix D2 (4. Ordnung)%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2i0                = -1;                                         %
f2i0                = -1*ones(numel(f2i0),n);                     %  
f2i                 = -30;                                        %
f2i                 = -30*ones(numel(f2i),n+2);                   %      
f2i1                = 16;                                         %
f2i1                = 16*ones(numel(f2i1),n+1);                   %
f2i2                = -1;                                         %
f2i2                = -1*ones(numel(f2i2),n);                     %
f2i3                = 16;                                         %
f2i3                = 16*ones(numel(f2i3),n+1);                   %
                                                                  %
D2_1                = (diag(f2i0,-2)+diag(f2i1,-1)+diag(f2i,0)+...%
                      diag(f2i3,1)+diag(f2i2,2));                 %

D2_1(1,1)           = 0;                                          %
D2_1(1,2)           = 0;                                          %
D2_1(1,3)           = 214;                                        %
D2_1(1,4)           = -156;                                       %
D2_1(1,5)           = +61;                                        %
D2_1(1,6)           = -10;                                        %

D2_1(2,1)           = 0;                                          %
D2_1(2,2)           = 0;                                          %
D2_1(2,3)           = 6;                                          %
D2_1(2,4)           = 4;                                          %
D2_1(2,5)           = -1;                                         %

D2_1(3,1)           = 0;                                          %
D2_1(3,2)           = 0;                                          % 
D2_1(4,2)           = 0;                                          % 

D2_1(n+2,n+2)         = 0;       
D2_1(n+2,n+1)         = 0;                                        %
D2_1(n+2,n)           = 214;                                      %
D2_1(n+2,n-1)         = -156;                                     %
D2_1(n+2,n-2)         = 61;   
D2_1(n+2,n-3)         = -10; 

D2_1(n+1,n+2)         = 0;  
D2_1(n+1,n+1)         = 0;  
D2_1(n+1,n)           = 6;                                        %
D2_1(n+1,n-1)         = 4;                                        %
D2_1(n+1,n-2)         = -1;   

D2_1(n,n+2)           = 0;                                        %
D2_1(n,n+1)           = 0;                                        %
D2_1(n-1,n+1)          = 0;                                       %
                                                                  %
D2                  = 1/(12*delta_y^2)*D2_1;                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Berechnung der Ableitungsmatrix D4 (4. Ordnung)%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4i0                = -2;                                         %
f4i0                = -2*ones(numel(f4i0),n-1);                   %
f4i                 = 112;                                        %
f4i                 = 112*ones(numel(f4i),n+2);                   %
f4i1                = -78;                                        %
f4i1                = -78*ones(numel(f4i1),n+1);                  %  
f4i2                = 24;                                         %
f4i2                = 24*ones(numel(f4i2),n);                     %
f4i4                = -78;                                        %
f4i4                = -78*ones(numel(f4i4),n+1);                  %
f4i5                = 24;                                         %
f4i5                = 24*ones(numel(f4i5),n);                     %
f4i6                = -2;                                         %
f4i6                = -2*ones(numel(f4i6),n-1);                   %
                                                                  %
D4_1                = diag(f4i0,-3)+diag(f4i2,-2)+...             %
                      diag(f4i1,-1)+diag(f4i,0)+diag(f4i4,+1)+... %
                      diag(f4i5,+2)+diag(f4i6,+3);                %
                                                                  %
D4_1(1,1)           = 0;                                          %
D4_1(1,2)           = 0;                                          %
D4_1(1,3)           = 1704;                                       %
D4_1(1,4)           = -2438;                                      %
D4_1(1,5)           = 2112;                                       %
D4_1(1,6)           = -1110;                                      %
D4_1(1,7)           = 328;                                        %
D4_1(1,8)           = -42;                                        %

D4_1(2,1)           = 0;                                          %
D4_1(2,2)           = 0;                                          %
D4_1(2,3)           = 342;                                        %
D4_1(2,4)           = -368;                                       %
D4_1(2,5)           = 222;                                        %
D4_1(2,6)           = -72;                                        %
D4_1(2,7)           = 10;                                         %

D4_1(3,1)           = 0;                                          %
D4_1(3,2)           = 0;                                          %
D4_1(3,3)           = 42;                                         %
D4_1(3,4)           = -8;                                         %
D4_1(3,5)           = -18;                                        %
D4_1(3,6)           = 12;                                         %
D4_1(3,7)           = -2;                                         %

D4_1(4,1)           = 0;                                          %
D4_1(4,2)           = 0;                                          %

D4_1(5,2)           = 0;                                          %

D4_1(n+2,n+2)         = 0;                                        %
D4_1(n+2,n+1)         = 0;                                        %
D4_1(n+2,n)           = 1704;                                     %
D4_1(n+2,n-1)         = -2438;                                    %
D4_1(n+2,n-2)         = 2112;                                     %
D4_1(n+2,n-3)         = -1110;                                    %
D4_1(n+2,n-4)         = 328;                                      %
D4_1(n+2,n-5)         = -42;                                      %

D4_1(n+1,n+2)           = 0;                                      %
D4_1(n+1,n+1)           = 0;                                      %
D4_1(n+1,n)            = 342;                                     %
D4_1(n+1,n-1)          = -368;                                    %
D4_1(n+1,n-2)          = 222;                                     %
D4_1(n+1,n-3)          = -72;                                     %
D4_1(n+1,n-4)          = 10;                                      %
 
D4_1(n,n+2)         = 0;                                          %
D4_1(n,n+1)         = 0;                                          %
D4_1(n,n)           = 42;                                         %
D4_1(n,n-1)         = -8;                                         %
D4_1(n,n-2)         = -18;                                        %
D4_1(n,n-3)         = 12;                                         %
D4_1(n,n-4)         = -2;                                         %

D4_1(n-1,n+2)         = 0;                                        %
D4_1(n-1,n+1)         = 0;                                        %

D4_1(n-2,n+1)         = 0;                                        %

D4                  = 1/(12*delta_y^4)*D4_1;                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U0_norm         = zeros(1,n); % Definition der dimensionslosen
                              % Grundstr�mung U0_norm(y_norm)

U0_norm2        = zeros(1,n); % Definition der zweiten Ableitung
                              % nach y_norm der dimensionslosen
                              % Grundstr�mung: U0_norm2(y_norm)
                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung des Blasius Profils

if blasius_flag==1
    options = odeset('AbsTol',1e-15,'RelTol',1e-12);
    [eta,f]=ode45(@BlasiusFunc, [0,100],[0 0 0.332],options);%u'(0)=0.3320430
    
    U=f(:,2);
    Y=eta/5.0*delta_delta1;
    U0_norm(1:n+2)=spline(Y,U,y_norm(1:n+2)); %Interpol. der Werte auf y_norm
else
    U0_norm(1:n+2)=spline(DNS_Y,DNS_U,y_norm(1:n+2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Berechnen der zweiten Ableitung nach y_norm der dimensionslosen %
% Grundstr�mung U0_norm (4. Ordnung)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=3:1:(n)                                                  
    U0_norm2(i)     = (-1*U0_norm(i-2)+16*U0_norm(i-1)-...
        30*U0_norm(i)+16*U0_norm(i+1)-1*U0_norm(i+2))/(12*delta_y^2);
end

U0_norm2(1)     = (45*U0_norm(1)-154*U0_norm(2)+214*U0_norm(3)-...
    156*U0_norm(4)+61*U0_norm(5)-10*U0_norm(6))/(12*delta_y^2);
U0_norm2(2)     = (11*U0_norm(1)-20*U0_norm(2)+6*U0_norm(3)+...
    4*U0_norm(4)-1*U0_norm(5))/(12*delta_y^2);
U0_norm2(n+1)   = (-1*U0_norm(n-2)+4*U0_norm(n-1)+6*U0_norm(n)-...
    20*U0_norm(n+1)+11*U0_norm(n+2))/(12*delta_y^2);
U0_norm2(n+2)     = (-10*U0_norm(n-3)+61*U0_norm(n-2)-...
    156*U0_norm(n-1)+214*U0_norm(n)-154*U0_norm(n+1)+...
    45*U0_norm(n+2))/(12*delta_y^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                              
                                  
fprintf('alpha-Intervall: %g -> %g\n',alpha_start,alpha_end)  ;                                

for alpha=alpha_start:alpha_schritt:alpha_end % Schleife f�r den Durchlauf 
                                              % mehrerer Wellenzahlen alpha

    fprintf('Untersuche Re_delta1=%f, alpha=%f - (%d/%d)\n',...
        Re_delta1,alpha,l,anzahl_alpha_schritte);
        
    % Berechnen der Vektoren der Koeffizienten a0 und a2 f�r die %%%%%%
    % Koeffizientenmatrizen A0 und A2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:1:n+2                                                     %
                                                                      %
        a0(i)           = -(alpha^3)*U0_norm(i)-alpha*U0_norm2(i)+... %
                          1i*(alpha^4/Re_delta1);                     %
        a2(i)           = alpha*U0_norm(i)-1i*(2*alpha^2/Re_delta1);  %
                                                                      %
    end                                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Besetzen der Koeffizientenmatrix A0 mit den berechneten      %%%%
    % Koeffizienten a0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposA0              = 1:n+2;                                      %
    jposA0              = 1:n+2;                                      %
    A0K                 = a0(1,:);                                    %
    A0                  = sparse (iposA0,jposA0,A0K,n+2,n+2);         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Besetzen der Koeffizientenmatrixatrix A2 mit den berechneten %%%%
    % Koeffizienten a0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposA2              = 1:n+2;                                      %
    jposA2              = 1:n+2;                                      %
    A2K                 = a2(1,:);                                    %
    A2                  = sparse (iposA2,jposA2,A2K,n+2,n+2);         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Erstellen der Koeffizientenmatrixatrix A4 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposA4              = 1:n+2;                                      %
    jposA4              = 1:n+2;                                      %
    A4K                 = 1i*(1/Re_delta1);                           %
    A4                  = sparse (iposA4,jposA4,A4K,n+2,n+2);         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Erstellen der Koeffizientenmatrixatrix B0 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposB0              = 1:n+2;                                      %
    jposB0              = 1:n+2;                                      %
    B0K                 = -(alpha^2);                                 %
    B0                  = sparse (iposB0,jposB0,B0K,n+2,n+2);         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Erstellen der Koeffizientenmatrixatrix B2 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposB2              = 1:n+2;                                      %
    jposB2              = 1:n+2;                                      %
    B2K                 = 1;                                          %
    B2                  = sparse (iposB2,jposB2,B2K,n+2,n+2);         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Berechnen der zusammengefassten Matrizen A und B zum Erreichen %%
    % der allgemeinen Form eines generalisierten Eigenwertproblems %%%%
    A                   = A0+A2*D2+A4*D4;                             %
    B                   = B0+B2*D2;                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [~,omega]         = eig(A,B); % L�sen des generalisierten
                                    % Eigenwertproblems mit dem
                                    % Eigenvektor phi und dem Eigenwert
                                    % omega

    im_omega            = diag(imag(omega)); % Auslesen der 
                                             % Imagin�rteile der
                                             % komplexen St�rfrequenz
                                             % omega: im_omega

    re_omega            = diag(real(omega)); % Auslesen der 
                                             % Realteile der komplexen
                                             % St�rfrequenz omega:
                                             % re_omega

    % Auslesen der Werte fuer Re_delta1, alpha und re_omega bei denen %%
    % die Anfachungsrate im_omega groeer als Null ist %%%%%%%%%%%%%%%%%
    
    number_amplified_eigenvalues_for_alpha = 0;
    
    for i=1:1:n+2                                                     %
                                                                      %
        if im_omega(i) > 0 && re_omega(i) > 0                         %
            number_amplified_eigenvalues=number_amplified_eigenvalues+1;
            number_amplified_eigenvalues_for_alpha = number_amplified_eigenvalues_for_alpha+1;
            inst_re_omega(number_amplified_eigenvalues,number_amplified_eigenvalues_for_alpha)  = re_omega(i);                        %
            inst_im_omega(number_amplified_eigenvalues,number_amplified_eigenvalues_for_alpha)  = im_omega(i);
            inst_alpha(number_amplified_eigenvalues,number_amplified_eigenvalues_for_alpha)     = alpha;                              %
            fprintf('angefachter Eigenwert gefunden: %f + i*%f\n', ...
                re_omega(i), im_omega(i));
            
            %Diese Abfrage prueft, ob der Eigenwert richtig ist (empirisch)
            %if (abs(re_omega(i) - alpha)/alpha>(0.001/alpha))
       
             omega_wielandt_alt_vektor(number_amplified_eigenvalues)=...
                    re_omega(i)+1i*im_omega(i);                       %
            %end
        end                                                           %
                                                                      %
    end
                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prüfen, ob mehrere Eigenwerte gefunden wurden und den am nächsten zu
% omega_alt liegenden auswählen

    if number_amplified_eigenvalues_for_alpha > 1
        for j=2:number_amplified_eigenvalues_for_alpha
            if inst_re_omega(number_amplified_eigenvalues,1)-omega_alt > inst_re_omega(number_amplified_eigenvalues,j)-omega_alt
                inst_re_omega(number_amplified_eigenvalues,1)=inst_re_omega(number_amplified_eigenvalues,j);
                inst_alpha(number_amplified_eigenvalues,1)=inst_alpha(number_amplified_eigenvalues,j);
                inst_im_omega(number_amplified_eigenvalues,1)=inst_im_omega(number_amplified_eigenvalues,j);
            end
        end
        omega_alt=inst_re_omega(number_amplified_eigenvalues,1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Aufruf der raeumlichen Loesung (wielandt) %%%%%%%%%%%%%%%%%%%%%
    if spatial_flag==1
       
       if number_amplified_eigenvalues > 1
           fprintf('ACHTUNG !!! Mehr als ein angefachter Eigenwert!\n');
       end
       
       for eigenvalue_iter=1:1:number_amplified_eigenvalues
           
           omega_wielandt_alt=omega_wielandt_alt_vektor(eigenvalue_iter);
           
           if imag(omega_wielandt_alt)>1.0e-5
               fprintf('alpha: %d/%d | angefachter Eigenwert :%d/%d\n',l,...
                   iterations_alpha,eigenvalue_iter,number_amplified_eigenvalues);
               
               [omega_wielandt_neu,alpha_wielandt_neu]=LST_NACA_spatial(alpha,omega_wielandt_alt, U0_norm, U0_norm2, D2, D4, Re_delta1) ;
               
               spatial_omega_r(k2)=real(omega_wielandt_neu);
               spatial_alpha_i(k2)=imag(alpha_wielandt_neu);
               spatial_alpha_r(k2)=real(alpha_wielandt_neu);
               k2=k2+1;
               
           end
       end
    end  
%%%%%%%%%%%%%%%%%raeumliche Loesung Ende%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

       l = l+1; % Erhoehen des Zaehlers fuer die Wellenzahl
       
end    %Schleife alpha Ende




inst_re_omega_final        = zeros(1,number_amplified_eigenvalues); % Definition des Vektors f�r die
                                         % Kreisfrequenz der 
                                         % St�rungspartialwelle der
                                         % ausschlie�lich Werte enth�lt bei
                                         % denen Instabilit�ten auftreten
                                         
inst_alpha_final           = zeros(1,number_amplified_eigenvalues); % Definition des Vektors f�r die
                                         % Reynolds-Zahl der ausschlie�lich
                                         % Werte enth�lt bei denen
                                         % Instabilit�ten auftreten
                                                                                                        
inst_c                      = zeros(1,number_amplified_eigenvalues); % Definition des Vektors f�r die 
                                          % Phasengeschwindigkeit der 
                                          % Partialwelle bei der 
                                          % Instabilit�ten auftreten:
                                          % inst_c
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%Gaster Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if number_amplified_eigenvalues > 0
    
    for i=2:number_amplified_eigenvalues-1
        dw_da(i)=((-1/2)*inst_re_omega(i-1,1)+(1/2)*inst_re_omega(i+1,1))/(inst_alpha(i+1,1)-inst_alpha(i,1));
    end
    dw_da(1)=((-3/2)*inst_re_omega(1,1)+(2)*inst_re_omega(2,1)+(-1/2)*inst_re_omega(3,1))/(inst_alpha(i+1,1)-inst_alpha(i,1));
    dw_da(number_amplified_eigenvalues)=dw_da(number_amplified_eigenvalues-1);
    
    alpha_spatial_gaster_i=inst_im_omega(:,1).*(-1)./transpose(dw_da);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% Berechnen der Phasengeschwindigkeit der Partialwelle bei der %%%%%%%%%
% Instabilit�ten auftreten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for x=1:1:number_amplified_eigenvalues                             %
                                                                       %
        inst_c(x) = inst_re_omega(x,1)/inst_alpha(x,1);                %
                                                                       %
    end                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    inst_F                      = zeros(1,number_amplified_eigenvalues); 
                                           % Definition des Vektors f�r die
                                           % dimensionslose Frequenz inst_F
                                           % bei der Instabilit�ten
                                           % auftreten
    
    d                           = 1;       % Initialisierung der Laufkoordinate
                                           % f�r die dimensionslose Frequenz inst_F
    
    
% Berechnen der dimensionslosen Frequenz inst_F bei der Instabilit�ten %
% auftreten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for x=1:1:number_amplified_eigenvalues                                                            %
        %
        inst_F(d) = (inst_alpha(x,1)/Re_delta1)*...      %
            inst_c(x)*(10^6);                                      %
        d         = d+1;                                                   %
        %
    end                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Plotten der zeitlichen Lösung
    
    %     figure,plot( re_omega,im_omega,'*');
    %
    %     axis([0,0.45,-.1,0.05]);
    %     xlabel('w_{r}');                                              %
    %     ylabel('w_{i}');                                               %
    %     title('eigenwerte');
    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Plotten von alpha_i (nach Gaster) über Omega_re
    
    f_real=inst_re_omega(:,1)/2/pi*(U_ref_DNS/L_ref);
    alpha_i = (-1)*alpha_spatial_gaster_i/L_ref;
    
    % figure,plot(f_real,alpha_i,'*');
    %
    % %axis([0,0.45,-.1,0.05]);
    % xlabel('Frequenz [Hz]');                                              %
    % ylabel('-alpha_i (gaster)');                                               %
    % title('alpha_i nach Gaster');
    
else
    alpha_i=0;
    f_real=0;
    fprintf('keinen instabilen Eigenwert gefunden!\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Plotten der Geschwindigkeitsverteilung U0_norm in Abh�ngikeit der %%%%%
 % dimensionslosen Koordinate y_norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure,plot(U0_norm(1:(n-1)),y_norm(2:n),'*');                                     
% xlabel('U/u_{\infty}'); ylabel('y/\delta_1');                           %
% title('Geschwindigkeitsverteilung-DNS');                                    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Plotten der zweiten Ableitung der Geschwindigkeitsverteilung U0_norm %%
% % in Abh�ngikeit der dimensionslosen Koordinate y_norm %%%%%%%%%%%%%%%%%%
% plot(U0_norm2(2:n),y_norm(2:n));                                        %
% xlabel('\partial^2{U}/\partial{y}^2'); ylabel('y/\delta_1');            %
% title('2. Ableitung der Geschwindigkeitsverteilung');                   %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Raeumliche Anregung
% Plotten der imaginären Wellenzahl in Abh�ngigkeit von der 
% realen Kreisfrequenz %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if spatial_flag==1
    figure;
    plot(spatial_omega_r,spatial_alpha_i, '--rs','LineWidth',2,...
                   'MarkerEdgeColor','k',...
                   'MarkerFaceColor','g',...
                   'MarkerSize',10)

    xlabel('\omega_r');                                              %
    ylabel('\alpha_i');                                               %
    title('Raeumliche Anregung'); 
    set(gca,'YDir','reverse');
end

% Plotten der imaginären Wellenzahl in Abh�ngigkeit von der 
% realen Wellenzahl %%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if spatial_flag==1
    figure;
    plot(spatial_alpha_r,spatial_alpha_i ,'--rs','LineWidth',2,...
                   'MarkerEdgeColor','k',...
                   'MarkerFaceColor','g',...
                   'MarkerSize',10)

    xlabel('\alpha_r');                                              %
    ylabel('\alpha_i');                                               %
    title('Raeumliche Anregung'); 
    set(gca,'YDir','reverse');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frequenz=spatial_omega_r/(2*pi)*U_inf/L_ref;
%frequenz=spatial_omega_r*U_inf^2/(2*pi*nu*Re_delta1);

% frequenz=spatial_omega_r/(2*pi)*U_inf/(L_ref*DNS_delta1);
% spatial_alpha_i_dim=spatial_alpha_i/(L_ref*DNS_delta1);
% figure,plot(frequenz,spatial_alpha_i_dim, '--rs','LineWidth',2,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',10)
%                      %
% xlabel('Frequenz[Hz]');                                              %
% ylabel('\alpha_i[rad/m]');                                               %
% title('Raeumliche Anregung'); 
% set(gca,'YDir','reverse');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Plotten der Phasengeschwindigkeit in Abh�ngigkeit von der %%%%%%%%%%%%%
% % Reynolds-Zahl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(inst_Re_delta1_final,inst_c,'*');                                  %
% xlabel('Re_{{\delta}_1}');                                              %
% ylabel('c/u_{\infty}');                                                 %
% title('Indifferenzkurve');                                              %
% set(gca, 'xscale', 'log');                                              %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Plotten der dimensionslosen Frequenz inst_F in Abh�ngigkeit von der %%%
% % Reynolds-Zahl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(inst_Re_delta1_final,inst_F,'*');                                  %
% xlabel('Re_{{\delta}_1}');                                              %
% ylabel('\omega_r\nu/u_{\infty}^2*10^6');                                %
% title('Indifferenzkurve');                                              %
% set(gca, 'xscale', 'log');                                              %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plotten der dimensionslosen Frequenz inst_F in Abh�ngigkeit von der %%%
% % Reynolds-Zahl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(inst_Re_delta1_final,inst_F,'*');                                  %
% xlabel('Re_{{\delta}_1}');                                              %
% ylabel('\omega_r\nu/u_{\infty}^2*10^6');                                %
% title('Indifferenzkurve');                                              %
% set(gca, 'xscale', 'log');                                              %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
