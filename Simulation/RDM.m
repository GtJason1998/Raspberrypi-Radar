%% Initial operation
close all;
clc;
tarR = [50 90];    %target range
tarV = [3 20];     %target velocity
c = 3*10^8;
f0 = 24.05*10^9;
T = 100*10^(-6);   %chirp Sweep Time
B = 400*10^6;
L = 600;            %slow-time dimension,num of chirps
N = 600;           %fast-time dimension,num of samples
%% generate receive signal
S1 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S1(l,n) = exp(1i*2*pi*(((2*B*(tarR(1)+tarV(1)*T*l)/(c*T)+(2*f0*tarV(1))/c)*T/N*n+((2*f0)*(tarR(1)+tarV(1)*T*l))/c)));
    end
end

S2 = zeros(L,N);
for l = 1:L
    for n = 1:N
        S2(l,n) = exp(1i*2*pi*(((2*B*(tarR(2)+tarV(2)*T*l)/(c*T)+(2*f0*tarV(2))/c)*T/N*n+((2*f0)*(tarR(2)+tarV(2)*T*l))/c)));
    end
end

sigReceive = S1+S2;
%% range fft processing
hanning1 = hanning(N,'periodic');
sigRWin = zeros(L,N);
for ii = 1:L
    sigRWin(ii,:) = hanning1'.*sigReceive(ii,:);
end
sigRfft = zeros(L,N);
for ii = 1:L
    sigRfft(ii,:) = fft(sigRWin(ii,:),N);
end
%% doppler fft processing
hanning2 = hanning(L,'periodic');
sigDWin = zeros(L,N);
for ii = 1:N
    sigDWin(:,ii) = hanning2.*sigRfft(:,ii);
end
sigDfft = zeros(L,N);
for ii = 1:N
    sigDfft(:,ii) = fft(sigDWin(:,ii),L);
    %sigDfft(:,ii) = fft(sigReceive(:,ii),NumDFFT);
end
%% Visualization
figure,image(abs(sigRfft)),title('Range - FFT')
figure,mesh(abs(sigRfft)),title('Range-FFT');
figure,image(abs(sigDfft)),title('RDM');xlabel('Range-FFT'),ylabel('Doppler-FFT');
figure,mesh(abs(sigDfft)),title('Range Dopple - FFT');