%Testsignal erstellen

Fs = 1000;                       % Sampling frequency
T = 1/Fs;                     % Sample time, Länge eines Zeitschritts
L = 1000;                     % Length of signal 
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = (3e-5)*sin(2*pi*7*t) + (1e-5)*sin(2*pi*20*t); 
y = x ;%+ 5*randn(size(t));     % Sinusoids plus noise

%Datei schreiben

fileID = fopen('/home/christiangscheidle/Dokumente/test/PressureHistory_Points_All.dat','w');

fprintf(fileID,' TITLE = "testsignal"\n');
fprintf(fileID,'VARIABLES = "t [s]" "PH_0_x-1_y-2_z0"\n');
fprintf(fileID,sprintf('ZONE T="testsignal", F=BLOCK, I=%d, DT=(SINGLE)\n',L));
for i=1:length(t)
    fprintf(fileID,sprintf('%d\n',t(i)));
end
for i=1:length(t)
    fprintf(fileID,sprintf('%d\n',y(i)));
end

fclose(fileID);
%Signal plotten

figure;
plot(Fs*t(1:end),y(1:end));
title('Signal Corrupted with Zero-Mean Random Noise');
xlabel('sample');
