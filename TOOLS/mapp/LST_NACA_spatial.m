function [omega_wielandt_neu,alpha_wielandt_neu]=LST_NACA_spatial(alpha,omega_wielandt_alt,U0_norm, U0_norm2, D2, D4, Re_delta1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raeumliche Loesung %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n=length(U0_norm);
alpha_wielandt_alt=alpha;
alpha_wielandt_neu=alpha_wielandt_alt+1i*10e-8;
TolMax = 10^-8 ;      %eps
TolOmegai = 10^-10;   %imag_Omega -> 0, eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a0                  = zeros(1,n); % Definition des Vektors der
% Koeffizienten a0 f�r die
% Koeffizientenmatrix A0

a2                  = zeros(1,n); % Definition des Vektors der
% Koeffizienten a2 f�r die
% Koeffizientenmatrix A2

% Besetzen der Koeffizientenmatrixatrix A0 mit den berechneten %%%%
% Koeffizienten a0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iposA0              = 1:n;                                        %
jposA0              = 1:n;                                        %
A0K                 = a0(1,:);                                    %
A0                  = sparse (iposA0,jposA0,A0K,n,n);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Besetzen der Koeffizientenmatrixatrix A2 mit den berechneten %%%%
% Koeffizienten a0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iposA2              = 1:n;                                        %
jposA2              = 1:n;                                        %
A2K                 = a2(1,:);                                    %
A2                  = sparse (iposA2,jposA2,A2K,n,n);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Erstellen der Koeffizientenmatrixatrix A4 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iposA4              = 1:n;                                        %
jposA4              = 1:n;                                        %
A4K                 = 1i*(1/Re_delta1);                           %
A4                  = sparse (iposA4,jposA4,A4K,n,n);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Erstellen der Koeffizientenmatrixatrix B2 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iposB2              = 1:n;                                        %
jposB2              = 1:n;                                        %
B2K                 = 1;                                          %
B2                  = sparse (iposB2,jposB2,B2K,n,n);             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Newton-Verfahren
break_newton=0;
for N=1:1:100
    % Berechnen der Vektoren der Koeffizienten a0 und a2 f�r die %%%%%%
    % Koeffizientenmatrizen A0 und A2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:1:n                                                       %
        %
        a0(i) = -(alpha_wielandt_neu^3)*U0_norm(i)-...                %
            alpha_wielandt_neu*U0_norm2(i)+...                        %
            1i*(alpha_wielandt_neu^4/Re_delta1);                      %
        a2(i) = alpha_wielandt_neu*U0_norm(i)-1i*...                  %
            (2*alpha_wielandt_neu^2/Re_delta1);                       %
        %
    end                                                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Erstellen der Koeffizientenmatrixatrix B0 %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iposB0              = 1:n;                                        %
    jposB0              = 1:n;                                        %
    B0K                 = -(alpha_wielandt_neu^2);                    %
    B0                  = sparse (iposB0,jposB0,B0K,n,n);             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A=(A0+A2*D2+A4*D4) / (B0+B2*D2);
    E=eye(length(A));
    x0= ones(length(A),1);      %Startwert 'Wielandt [1;1;1]
    Bw = (A-omega_wielandt_alt*E);
    
    break_wielandt=0;
    % Wielandt
    for W=1:1:100
        x1 = Bw\x0;
        [zeile,spalte] = find(abs(x1)==max(abs(x1)));
        if length(zeile)>1
            zeile=zeile(1);
            spalte=spalte(1);
        end
        Max1 = x1(zeile,spalte);
        x1Norm = Max1^-1 * x1;
        
        if W > 1
            if abs(Max1)-abs(Max0) < TolMax
                break_wielandt=1;
                break % Schleife Wielandt
            end
        end
        
        x0 = x1Norm;
        Max0 =Max1 ;
    end % Schleife Wielandt
    %% ------------------------
    omega_wielandt_neu = omega_wielandt_alt + 1/Max1;
    
    %% Newton-Verfahren
    if abs(imag(omega_wielandt_neu))< TolOmegai
        break_newton=1;
        break % Schleife Newton
    else
        grad = (imag(omega_wielandt_neu)-imag(omega_wielandt_alt))/...
            (imag(alpha_wielandt_neu)-imag(alpha_wielandt_alt));
        alpha_wielandt_alt = alpha_wielandt_neu ;
        alpha_wielandt_neu = alpha_wielandt_alt - 1i * ...
            imag(omega_wielandt_neu)/grad;
        omega_wielandt_alt = omega_wielandt_neu;
    end
end % Schleife Newton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%msg = sprintf('Berechne für Re_delta1: %d/%d und alpha: %d/%d |...
%angefachte Eigenwerte :%d',s,iterations_Re_delta1,l,...
%iterations_alpha,number_amplified_eigenvalues);
%fprintf([reverseStr, msg]);
%reverseStr = repmat(sprintf('\b'), 1, length(msg));


if real(omega_wielandt_neu) > 0.0000
    %if break_newton==1 && break_wielandt==1 && ...
    %imag(alpha_wielandt_neu)<0.0
    if imag(alpha_wielandt_neu)<0.0
        fprintf('\nSpatial: ')
        fprintf('omega = %f %f*i | ',real(omega_wielandt_neu), ...
            imag(omega_wielandt_neu))
        fprintf('alpha= %f  %f*i \n\n',real(alpha_wielandt_neu),...
            imag(alpha_wielandt_neu))
    end
end


end
