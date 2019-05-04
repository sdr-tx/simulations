clear, clc, close all;
%% Variables
fc = 10e6;       Tc = 1/fc;
fpwm = 1e6;     Tpwm = 1/fpwm;

f_clk = 100e7;    T_clk = 1/f_clk;

%%
Ncyc_pwm = 5;
Ncyc_cry = fc/fpwm* Ncyc_pwm;
Ncyc_clk = f_clk/fc * Ncyc_cry;

%%
[t1, carrier] = cuadrada(Ncyc_cry, 0.5, fc, f_clk);
[pwm_1] = cuadrada_con_desfasaje(Ncyc_pwm, 1, 0, 0.75, fpwm, f_clk);
sig1 = pwm_1.*carrier;

pwm2 = cuadrada_con_desfasaje(Ncyc_pwm, 1, -.25, 0.25, fpwm, f_clk);
sig2 = pwm2.*carrier;

pwm3 = cuadrada_con_desfasaje(Ncyc_pwm, 1, .25, 0.5, fpwm, f_clk);
sig3 = pwm3.*carrier;

sig = [zeros(1,2000) sig1 sig2 sig3];
t = 0:T_clk:(length(sig)-1)*T_clk;

Num = load('coeficientes.fcf');
sig_filt = filter(Num,1,sig);

%%
Modulation_spect = abs(fft(sig)) / length(sig);
f_modu = f_clk/2 * linspace(0, 1, length(Modulation_spect)/2+1);

Modulation_spect_filt = abs(fft(sig_filt)) / length(sig_filt);

figure
h1=subplot(211),plot(f_modu, 2*Modulation_spect(1:length(Modulation_spect)/2+1)); title('Signal Spectrum')
h2=subplot(212),plot(f_modu, 2*Modulation_spect_filt(1:length(Modulation_spect_filt)/2+1)); 
title('Filtered Signal Spectrum')
linkaxes([h1,h2],'x');
%%
z = qamdemod(sig_filt,4);
figure,
h3 = subplot(311),plot(t,sig), ylim([-0.5,1.5]), grid
title('PWMed QAM: Symb1(D=0.75, Ph=0); Symb2 (D=0.25, Ph=-0.25)')

h4 = subplot(312),plot(t,sig_filt), ylim([-1.5,1.5]), grid
title('PWMed QAM filtered')

h5 = subplot(313),plot(t, z),  ylim([-0.5,3]),title('Demodulación')
linkaxes([h3,h4,h5],'x')
scatterplot(z)

