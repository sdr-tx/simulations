clc, , close all
%%
% genero 5 ciclos por s�mbolo
fase1 = 0;
fase2 = 0.25;
fase3 = 0.5;
fase4 = -0.25;

fpwm = 1e3; Tpwm = 1/fpwm;
fs = 1e6;   dt=1/fs;

ncyc = 10;

t = 0: dt: 4*ncyc*Tpwm-dt;

duty = 0.5;


pwm1 = cuadrada_con_desfasaje(ncyc,1,fase1,0.6,fpwm,fs);
pwm2 = cuadrada_con_desfasaje(ncyc,1,fase2,duty,fpwm,fs);
pwm3 = cuadrada_con_desfasaje(ncyc,1,fase3,duty,fpwm,fs);
pwm4 = cuadrada_con_desfasaje(ncyc,1,fase4,duty,fpwm,fs);

signal = [pwm1 pwm2 pwm3 pwm4];

%signal = [pwm1 pwm3];


figure
plot(t, signal); ylim([-0.5 1.5]); grid

%%
ventana = blackman(length(signal));

%signal_ventaneada = signal.*ventana';

Modulation_spect = abs(fft(signal)) / length(signal);
%Modulation_spect_ventaneado = abs(fft(signal_ventaneada)) / length(signal);

f_modu = fs/2 * linspace(0, 1, length(Modulation_spect)/2+1);
%%
num = load('coeficientes_QPSK_may_10.fcf');
filtered_signal = filter(num,1,signal);
figure
plot(t,filtered_signal);

%%

Modulation_spect_filtered = abs(fft(filtered_signal)) / length(signal);
%Modulation_spect_ventaneado = abs(fft(signal_ventaneada)) / length(signal);

f_modu = fs/2 * linspace(0, 1, length(Modulation_spect_filtered)/2+1);

figure
h1=subplot(211),plot(f_modu, 2*Modulation_spect(1:length(Modulation_spect)/2+1)); title('Signal Spectrum')
h2=subplot(212),plot(f_modu, 2*Modulation_spect_filtered(1:length(Modulation_spect_filtered)/2+1)); title('Signal Spectrum'), hold on
linkaxes([h1,h2],'x');

%%

I = filtered_signal .* cos(2*pi*fpwm*t);
Q = filtered_signal .* sin(2*pi*fpwm*t);
plot(t,I,t,Q)


%% DEMODULATION

i_filtered = filter(Num, 1, I);
q_filtered = filter(Num, 1, Q);

iq = (i_filtered) + 1i * (q_filtered);
% message = pskdemod(iq, 4);
message = qamdemod(iq, 8);


%% PLOT
figure(3);
subplot(121);   plot(t, i_filtered);        title("i in time");
subplot(122);   plot(t, q_filtered);        title("q in time");

scatterplot(iq);        title("i+q");
scatterplot(message);	title("message");

figure(4);
subplot(121);   plot(message, 'linewidth', 2);         title("message in time");
subplot(122);   histogram(message);    title("message histogram");

