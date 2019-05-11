%%  Modulaci�n QPSK:
%       - 4 se�ales de PWM con mismo duty (50%)
%       - el desfasaje de las mismas es de 0�, 90�, 180� y 270�
%       - frecuencia de muestreo de la se�al: fs = 250MHz (concordante con la de clock
%       de la FPGA)
%       - frecuencia de PWM: fpwm = 25MHz (la decima parte de fs)
%
%   Ensayos a realizar:
%       - Qu� cambia al variar la frecunecia de PWM?
%       - Qu� cambia al variar la cantidad de periodos de PWM enviados por
%       s�mbolo
%       - Ver c�mo reacciona el sistema con filtros con distintos anchos de
%       banda
%%
clc, , close all, clear
%%  Generaci�n de la se�ales en tiempo y filtrado

%Frecuencia de muestreo (250MHz)
fs = 250e6;     dt = 1/fs;

%Frecuencia de PWM (25MHz, fs/10)
fpwm = 25e6; Tpwm = 1/fpwm;

%Declaro las 4 fases distintas
fase1 = 0;
fase2 = -0.25;
fase3 = -0.5;
fase4 = -0.75;

%Cantidad de ciclos de cada PWM
ncyc = 100;

%Cantidad de s�mbolos a enviar
nsym = 8;

%Genero el eje de tiempo: 
% tfinal = cantidad de simbolos * ciclos de pwm * periodo del PWM
t = 0: dt: nsym*ncyc*Tpwm-dt;

%Duty de los PWM
duty = 0.5;

%Genero los nsym de PWM:
%La funci�n recibe CantidadDeCiclos, Amplitud, Fase, Duty, Frecuencia de
%PWM y Frecuencia de sampleo
pwm1 = cuadrada_con_desfasaje(ncyc,1,fase1,duty,fpwm,fs);
pwm2 = cuadrada_con_desfasaje(ncyc,1,fase2,duty,fpwm,fs);
pwm3 = cuadrada_con_desfasaje(ncyc,1,fase3,duty,fpwm,fs);
pwm4 = cuadrada_con_desfasaje(ncyc,1,fase4,duty,fpwm,fs);

%Se�al final: concatenaci�n de las PWM
signal = [pwm1 pwm3 pwm2 pwm4 pwm1 pwm2 pwm3 pwm4];

figure('name','Se�ales temporales')
plot(t, signal); ylim([-0.5 1.5]); grid; title('Se�al temporal')


%% Filtrado

%Cargo filtros
load('QPSK_filtro_FIR_fp24MHz_fc23MHz_at40dB.mat');
load('QPSK_filtro_IIR_fp24MHz_fc23MHz_at40dB.mat');

signal_filtered_FIR = filter(FIR_fc25MHz_fs23MHz,signal);
signal_filtered_IIR = filter(Hd,signal);

figure('name','Se�ales temporales')
h1=subplot(311),plot(t, signal); title('Se�al temporal'); ylim([-0.5 1.5]); grid;
h2=subplot(312),plot(t, signal_filtered_FIR); title('Se�al filtrada con FIR'); ylim([-1 1]); grid;
h3=subplot(313),plot(t, signal_filtered_IIR); title('Se�al filtrada con IIR'); ylim([-1 1]); grid;
linkaxes([h1,h2,h3],'x');


%% Espectros de las se�ales (filtrada y sin filtrar)

%Eje de frecuencia (medio eje)
f = fs/2 * linspace(0, 1, length(signal)/2+1);

%Expectro de la se�al filtrada y sin filtrar
signal_spect = abs(fft(signal)) / length(signal);
signal_spect_FIR = abs(fft(signal_filtered_FIR)) / length(signal);
signal_spect_IIR = abs(fft(signal_filtered_IIR)) / length(signal);

figure('name','Espectros')
h1=subplot(311), plot(f, 2*signal_spect(1:length(signal_spect)/2+1)); title('Signal Spectrum'),grid
h2=subplot(312), plot(f, 2*signal_spect_FIR(1:length(signal_spect)/2+1)); title('Signal Spectrum FIR'),grid
h3=subplot(313), plot(f, 2*signal_spect_IIR(1:length(signal_spect)/2+1)); title('Signal Spectrum IIR'),grid
linkaxes([h1,h2,h3],'x');

%% Descomposici�n en cuadratura

figure('name','se�ales en cuadratura')
I_FIR = signal_filtered_FIR .* cos(2*pi*fpwm*t);
Q_FIR = signal_filtered_FIR .* sin(2*pi*fpwm*t);
h1 = subplot(211), plot(t,I_FIR,t,Q_FIR), grid, title('I Q FIR sin LPF')

I_IIR = signal_filtered_IIR .* cos(2*pi*fpwm*t);
Q_IIR = signal_filtered_IIR .* sin(2*pi*fpwm*t);
h2 = subplot(212), plot(t,I_FIR,t,Q_FIR), grid, title('I Q IIR sin LPF')

linkaxes([h1 h2],'x')

%% DEMODULATION

%Cargo el filtro pasa bajos
load('QPSK_LPF_fp100kHz_fs25MHz');

%Filtro las se�ales I Q
I_FIR_filtered = filter(LPF_fp100kHz_fstop25MHz,I_FIR);
Q_FIR_filtered = filter(LPF_fp100kHz_fstop25MHz,Q_FIR);

I_IIR_filtered = filter(LPF_fp100kHz_fstop25MHz,I_IIR);
Q_IIR_filtered = filter(LPF_fp100kHz_fstop25MHz,Q_IIR);

figure
h1 = subplot(211); plot(t,I_FIR_filtered,t,Q_FIR_filtered), grid, title('I Q con FIR Filtrado')
h2 = subplot(212); plot(t,I_IIR_filtered,t,Q_IIR_filtered), grid, title('I Q con IIR Filtrado')

%% Genero la se�al en cuadratura que va al demodulador
iq_FIR = (I_FIR_filtered) + 1i * (Q_FIR_filtered);
iq_IIR = (I_IIR_filtered) + 1i * (Q_IIR_filtered);
message_FIR = pskdemod(iq_FIR, 4);
message_IIR = pskdemod(iq_IIR, 4);

scatterplot(iq_FIR);        title('I+Q FIR');
scatterplot(iq_IIR);        title('I+Q IIR');
%subplot(211), scatterplot(message_FIR);	title('Mensaje');
%subplot(212), scatterplot(message_IIR);	title('Mensaje');


figure;
subplot(121);   plot(t,message_FIR,'linewidth', 2); title('Mensaje en tiempo FIR');
subplot(122);   plot(t,message_FIR,'linewidth', 2); title('Mensaje en tiempo IIR');

subplot(121); hist(message_FIR);    title('Histograma FIR');
subplot(122); hist(message_FIR);    title('Histograma IIR');
