%%  Modulación QPSK:
%       - 4 señales de PWM con mismo duty (50%)
%       - el desfasaje de las mismas es de 0°, 90°, 180° y 270°
%       - frecuencia de muestreo de la señal: fs = 250MHz (concordante con la de clock
%       de la FPGA)
%       - frecuencia de PWM: fpwm = 25MHz (la decima parte de fs)
%
%   Ensayos a realizar:
%       - Qué cambia al variar la frecunecia de PWM?
%       - Qué cambia al variar la cantidad de periodos de PWM enviados por
%       símbolo
%       - Ver cómo reacciona el sistema con filtros con distintos anchos de
%       banda
%%
clc, , close all, clear
tic
%%  Generación de la señales en tiempo y filtrado

%Frecuencia de muestreo (250MHz)
fs = 2500e6;     dt = 1/fs;

%Frecuencia de PWM (25MHz, fs/10)
fpwm = 25e6; Tpwm = 1/fpwm;

%Declaro las 4 fases distintas
fase1 = 0;
fase2 = -0.25;
fase3 = -0.5;
fase4 = -0.75;

duty1 = 0.5;
duty2 = 0.25;

%Cantidad de ciclos de cada PWM
ncyc = 100;

%Cantidad de símbolos a enviar
nsym = 8;

%Genero el eje de tiempo: 
% tfinal = cantidad de simbolos * ciclos de pwm * periodo del PWM
t = 0: dt: nsym*ncyc*Tpwm-dt;

%Duty de los PWM
duty = 0.5;

%Genero los nsym de PWM:
%La función recibe CantidadDeCiclos, Amplitud, Fase, Duty, Frecuencia de
%PWM y Frecuencia de sampleo
pwm1 = cuadrada_con_desfasaje(ncyc,1,fase1,duty1,fpwm,fs);
pwm2 = cuadrada_con_desfasaje(ncyc,1,fase2,duty1,fpwm,fs);
pwm3 = cuadrada_con_desfasaje(ncyc,1,fase3,duty1,fpwm,fs);
pwm4 = cuadrada_con_desfasaje(ncyc,1,fase4,duty1,fpwm,fs);
pwm5 = cuadrada_con_desfasaje(ncyc,1,fase1,duty2,fpwm,fs);
pwm6 = cuadrada_con_desfasaje(ncyc,1,fase2,duty2,fpwm,fs);
pwm7 = cuadrada_con_desfasaje(ncyc,1,fase3,duty2,fpwm,fs);
pwm8 = cuadrada_con_desfasaje(ncyc,1,fase4,duty2,fpwm,fs);

%Señal final: concatenación de las PWM

signal = [pwm1 pwm5 pwm2 pwm6 pwm3 pwm7 pwm4 pwm8];

figure('name','Señales temporales')
plot(t, signal); ylim([-0.5 1.5]); grid; title('Señal temporal')



%% Filtrado

%Cargo filtros
load('QPSK_BPF_FIR_fp28MHz_fc50MHz_at60dB.mat');
load('QPSK_BPF_IIR_fp28MHz_fc50MHz_at60dB.mat');

signal_filtered_FIR = filter(BPF,signal);
signal_filtered_IIR = filter(BPF_IIR,signal);

figure('name','Señales temporales')
h1=subplot(311),plot(t, signal); title('Señal temporal'); ylim([-0.5 1.5]); grid;
h2=subplot(312),plot(t, signal_filtered_FIR); title('Señal filtrada con FIR'); ylim([-1 1]); grid;
h3=subplot(313),plot(t, signal_filtered_IIR); title('Señal filtrada con IIR'); ylim([-1 1]); grid;
linkaxes([h1,h2,h3],'x');


%% Espectros de las señales (filtrada y sin filtrar)

%Eje de frecuencia (medio eje)
f = fs/2 * linspace(0, 1, length(signal)/2+1);

%Expectro de la señal filtrada y sin filtrar
signal_spect = abs(fft(signal)) / length(signal);
signal_spect_FIR = abs(fft(signal_filtered_FIR)) / length(signal);
signal_spect_IIR = abs(fft(signal_filtered_IIR)) / length(signal);

figure('name','Espectros')
h1=subplot(311), plot(f, 2*signal_spect(1:length(signal_spect)/2+1)); title('Signal Spectrum'),grid
h2=subplot(312), plot(f, 2*signal_spect_FIR(1:length(signal_spect)/2+1)); title('Signal Spectrum FIR'),grid
[h, w] = freqz(BPF, length(t));
w = w(1:length(signal_spect)/2+1);
h3=subplot(313), plot(w * fs/2, abs(h(1:length(signal_spect)/2+1))); title('filter'),grid


linkaxes([h1,h2,h3],'x');

%% Descomposición en cuadratura

figure('name','señales en cuadratura')
I_FIR = signal_filtered_FIR .* cos(2*pi*fpwm*t);
Q_FIR = signal_filtered_FIR .* sin(2*pi*fpwm*t);
h1 = subplot(211), plot(t,I_FIR,t,Q_FIR), grid, title('I Q FIR sin LPF')

I_IIR = signal_filtered_IIR .* cos(2*pi*fpwm*t);
Q_IIR = signal_filtered_IIR .* sin(2*pi*fpwm*t);
h2 = subplot(212), plot(t,I_IIR,t,Q_IIR), grid, title('I Q IIR sin LPF')
linkaxes([h1 h2],'x')

%% DEMODULATION

%Cargo el filtro pasa bajos
load('QPSK_LPF_FIR_pf1MHz_fc5MHz_at60dB.mat');
%load('QPSK_IPF_FIR_pf1MHz_fc5MHz_at60dB.mat');


%Filtro las señales I Q
I_FIR_filtered = filter(LPF,I_FIR);
Q_FIR_filtered = filter(LPF,Q_FIR);

I_IIR_filtered = filter(LPF,I_IIR);
Q_IIR_filtered = filter(LPF,Q_IIR);

figure
h1 = subplot(211); plot(t,I_FIR_filtered,t,Q_FIR_filtered,'linewidth',2.5), grid, title('I Q con FIR Filtrado')
h2 = subplot(212); plot(t,I_IIR_filtered,t,Q_IIR_filtered,'linewidth',2.5), grid, title('I Q con IIR Filtrado')

%% Genero la señal en cuadratura que va al demodulador
iq_FIR = (I_FIR_filtered) + 1i * (Q_FIR_filtered);
iq_IIR = (I_IIR_filtered) + 1i * (Q_IIR_filtered);
message_FIR = qamdemod(iq_FIR, 16);
message_IIR = qamdemod(iq_IIR, 16);

scatterplot(iq_FIR);        title('I+Q FIR');
%scatterplot(iq_IIR);        title('I+Q IIR');
% subplot(211), scatterplot(message_FIR);	title('Mensaje');
% subplot(212), scatterplot(message_IIR);	title('Mensaje');


figure;
subplot(221);   plot(t,message_FIR,'linewidth', 2); title('Mensaje en tiempo FIR');
subplot(222);  plot(t,message_IIR,'linewidth', 2); title('Mensaje en tiempo IIR');
subplot(223); hist(message_FIR);    title('Histograma FIR');
subplot(224); hist(message_IIR);    title('Histograma IIR');
