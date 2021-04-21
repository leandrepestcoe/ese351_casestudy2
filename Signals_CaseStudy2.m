% Case Study 2
% Leandre Pestcoe and Julianne Wegmann
% ESE 351: Signals and Systems
% Created on: 4/20/21, Last Edited on: 4/21/21

%% Define Pulse Shape p(t)
Ts = 0.1; %symbol period (rate 1/Ts)
dt = 0.01; %sample period
w = 5*Ts; %width
t = -w:dt:w; %time vector
fs = 1/dt; %sample frequency

%sinc
sinc_p_t = sinc(t/Ts); %define sinc

%triangular pulse
triang_p_t = tripuls(t,w*2);

%decide which signal to use for simulation
p_t = sinc_p_t;

figure, plot(t,p_t), grid on;
xlabel('time (s)'), ylabel('p(t)'), title('Truncated Signal p(t)')

%% Nyquist Criteria

Nfft = 1024; %length of fft
f1 = (-fs/2:fs/Nfft:(fs/2));
f2 = (fs-2*pi/Ts:fs/Nfft:fs+(2*pi/Ts));
fft_signal = abs(fft(p_t,Nfft));
fft_signal1 = [fft_signal((length(fft_signal)/2)+1:end), fft_signal(1:(length(fft_signal)/2))];

figure
plot(f1(1:length(fft_signal1)),fft_signal1), grid on;
hold on
plot(f2(1:length(fft_signal1)),fft_signal1), grid on;
xlabel('frequency (Hz)'), ylabel('|Y(j\omega)|');

%% Noise-free PAM Signal y(t)
N = 10;
bits1 = 2*((rand(1,N)>0.5)-0.5);
bits2 = 2*((rand(1,N)>0.5)-0.5);
bits3 = 2*((rand(1,N)>0.5)-0.5);

x_t1 = zeros(1,N*(Ts/dt));
for i=1:length(bits1)
    x_t1((i-1)*(Ts/dt)+1)=bits1(i);
end
y_t1 = conv(x_t1,p_t);
t = (0:length(y_t1)-1)*dt;
x_new = zeros(1,length(y_t1));
x_new((length(p_t)+1)/2:(length(x_t1)+(length(p_t)+1)/2)-1) = x_t1;

figure();
plot(t,y_t1), grid on;
hold on
stem(t,x_new);
title('Transmitted Signal y(t)');
xlabel('Time[s]');
ylabel('y(t)');

%% Signal Modulation (Up-Conversion)
wc20 = 2*pi*20; %20 Hz modulation
mod_signal20 = y_t1.*cos(wc*t);

figure
subplot(2,1,1), plot(t,mod_signal), grid on;
xlabel('time (s)'), ylabel('y(t)cos(wc*t)'), title('Modulated Signal')

Nfft = 1024; %length of fft
f = (0:fs/Nfft:fs-fs/Nfft);
subplot(2,1,2), plot(f,abs(fft(mod_signal,Nfft))), grid on;
xlabel('frequency (Hz)'), ylabel('|Y(j\omega)|');

%% Noisy Recieved Signal r(t)
sigma = 0.5;
n_t = sigma*randn(1,length(mod_signal));
r_t = mod_signal + n_t;

plot(t,r_t(1:length(t)));
title('Received Signal r(t)');
xlabel('Time[s]');
ylabel('Amplitude');

%% Signal Demodulation (Down-Conversion)
demod_signal = r_t.*cos(wc*t);

figure
subplot(2,1,1), plot(t,demod_signal), grid on;
xlabel('time (s)'), ylabel('x_r(t)'), title('Demod (without LPF)')

f = (0:fs/Nfft:fs-fs/Nfft);
subplot(2,1,2), plot(f,abs(fft(demod_signal,Nfft))), grid on;
xlabel('frequency (Hz)'), ylabel('|X_r(j\omega)|')

%% Matched Filter Receiver

xn_tilda = ones(1,length(t));
p_neg = p_t(end:-1:1);
z_t = conv(r_t,p_neg);
for i=1:length(z_t)
    if z_t(i)<=0
        xn_tilda(i)=-1;
    end
end

%plot(t,xn_tilda(1:length(t)));
plot(z_t);

%% Noise Levels and Error Rates


