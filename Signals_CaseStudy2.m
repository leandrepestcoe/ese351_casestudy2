% Case Study 2
% Leandre Pestcoe and Julianne Wegmann
% ESE 351: Signals and Systems
% Created on: 4/20/21, Last Edited on: 4/20/21

%% Define pulse shape
Ts = 0.1; %symbol period (rate 1/Ts)
dt = 0.01; %sample period
w = 5*Ts; %width
t = -w:dt:w; %time vector
fs = 1/dt; %sample frequency

%sinc
%sinc = sinc(t/Ts); %define sinc

%triangular pulse
triang = tripuls(t,w*2);

%rectangular pulse
rect = zeros(1,length(t));
for i=(1:length(t)/2)
    rect(i)=1;
end

%decide which signal to use for simulation
signal = rect;

figure, plot(t,signal);
xlabel('time (s)'), ylabel('x(t)'), title('Truncated Signal')

%% Signal Modulation
wc = 2*pi*20; %20 Hz modulation
mod_signal = signal.*cos(wc*t);

figure
subplot(2,1,1), plot(t,mod_signal), grid on;
xlabel('time (s)'), ylabel('y(t)'), title('Modulated Signal')

Nfft = 1024; %length of fft
f = (0:fs/Nfft:fs-fs/Nfft);
subplot(2,1,2), plot(f,abs(fft(mod_signal,Nfft))), grid on;
xlabel('frequency (Hz)'), ylabel('|Y(j\omega)|');

%% Signal Demodulation
demod_signal = mod_signal.*cos(wc*t);

figure
subplot(2,1,1), plot(t,demod_signal), grid on;
xlabel('time (s)'), ylabel('x_r(t)'), title('Demod (without LPF)')

f = (0:fs/Nfft:fs-fs/Nfft);
subplot(2,1,2), plot(f,abs(fft(demod_signal,Nfft))), grid on;
xlabel('frequency (Hz)'), ylabel('|X_r(j\omega)|')

%%
% c) binary message 
N = 10;
x_n = 2*((rand(1,N)>0.5)-0.5);

% d) y(t)
t_p = 0.1;
time = -t_p:1/fs:((N-1)*2*t_p)+t_p;
x_t = zeros(1,length(time));
i = 1;
for j=11:10:length(time);
    x_t(j)=x_n(i);
    i = i+1;
    if i==N+1
        break 
    end
end
y_t = conv(x_t,p_t);
y_t_new = y_t(11:211);

% e) r(t) noise
sigma = 0.5;
n_t = sigma*randn(1,length(y_t_new));
r_t = y_t_new + n_t;

% f) Plots
% i) Pulse shape of p(t), and spectrum, P(w)
figure, plot(t,p_t), grid on, title('p(t) Pulse Shape'), xlabel('Time'), ylabel('Amplitude');
figure, plot(t,fft(p_t)), grid on, title('P(w)'), xlabel('Time'), ylabel('Amplitude');

% ii) Noise-free transmitted signal, y(t)
figure, plot(time,y_t_new), grid on, title('y(t), Noise-free Transmitted Signal'),...
    xlabel('Time'), ylabel('Amplitude');


% iii) Noisy received signal, r(t)
figure, plot(time,r_t), grid on, title('r(t), Noisy Received Signal'), xlabel('Time'), ylabel('Amplitude'), legend('Sigma=0.5, Bit Rate= 1/2*T_p');

% iv) Sent message, x_n, and decoded message, x_hat_n
% 1) Sign-based receiver 
i = 1;
for l=11:10:201
    if r_t(l) <= 0;
        x_n_received(i) = -1;
        i = i + 1;
    end
    if r_t(l) > 0;
        x_n_received(i) = 1;
        i = i + 1;
    end
    if i==N+1
        break 
    end
end
figure, hold on, grid on, plot(x_n,'--o'), plot(x_n_received), legend('x_n', 'x_n__received'), text(3.5,0.5,'Sigma=0.5, Bit Rate= 1/2*T_p, Error Rate = 0'),...
     title('Sign-based Receiver'),xlabel('bit number'),ylabel('bit value');

% error rate for sign-based receiver
correct = 0;
for k=1:10
    if x_n(k)==x_n_received(k)
        correct=correct+1;
    end
end
error_signbased = (10-correct)/10;

% 2) Matched Filter

% z_t
z_t = conv(r_t,p_t);
z_t_new = z_t(11:211);
i = 1;
for l=11:10:201
    if z_t_new(l) <= 0;
        x_n_received2(i) = -1;
        i = i + 1;
    end
    if z_t_new(l) > 0;
        x_n_received2(i) = 1;
        i = i + 1;
    end
    if i==N+1
        break 
    end
end
figure, hold on, grid on, plot(x_n,'--o'), plot(x_n_received2), legend('x_n', 'x_n__received2'), text(3.5,0.5,'Sigma=0.5, Bit Rate= 1/2*T_p, Error Rate = 0'),...
    title('Matched Filter'),xlabel('bit number'),ylabel('bit value');


% error rate for matched filter
correct = 0;
for k=1:10
    if x_n(k)==x_n_received2(k)
        correct=correct+1;
    end
end
error_matched = (10-correct)/10;