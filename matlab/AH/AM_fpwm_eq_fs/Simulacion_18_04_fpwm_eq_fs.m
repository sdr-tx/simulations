clear, clc, close all;
%% Variables
fc = 1e6;   Tc = 1/fc;
fs = 50e3;      Ts = 1/fs;

N = 8;
COMB = 2^N;

f_clk = fs*COMB;   T_clk = 1/f_clk;
f_pwm = fs;     T_pwm = 1/f_pwm;

%%
Ncyc = 3;
f_mod = 1e3;    w_mod = 2*pi*f_mod;     T_mod = 1/f_mod;

t = 0: Ts: Ncyc * T_mod-Ts;

signal = cos(w_mod*t);% + 2*cos(5*w_mod*t);

sig_Nbits = signal + abs(min(signal));
signal_2 = sig_Nbits / max(sig_Nbits) *(2^N-1);
sig_Nbits = round((sig_Nbits/max(sig_Nbits))*(2^N-1));

%% Generar PWM
pwm = 0;
for k = 1 : length(sig_Nbits)
    cant_1s = sig_Nbits(k) * f_clk/f_pwm / 2^N;
    cant_0s = 2^N - cant_1s;
    pwm = horzcat(pwm, ones(1, cant_1s), zeros(1, cant_0s));
end
pwm = pwm(2:end);
t_pwm = (0:length(pwm)-1)*T_clk;


%% PWM + PORTADORA

fc = f_clk/10;   Tc = 1/fc;
Ncyc_carrier = fc/f_mod * Ncyc;
[t_carrier, carrier] = cuadrada(Ncyc_carrier, 0.5, fc, f_clk);

Modulation = carrier.*pwm;



h1 = subplot(311); stem(t,sig_Nbits); title('Señal modulante');
h2 = subplot(312); plot(t_pwm,pwm); ylim([-0.5, 1.5]); title('PWM de la modulante')
h3 = subplot(313); plot(t_carrier, Modulation); ylim([-0.5, 1.5]), title('PWM + Carrier')
    
linkaxes([h1,h2,h3],'x');

%%

Modulation_spect = abs(fft(Modulation)) / length(Modulation);
f_modu = f_clk/2 * linspace(0, 1, length(Modulation_spect)/2+1);

coef = mi_filtro_coef(f_clk, fc-40*f_mod, fc-10*f_mod, fc+10*f_mod, fc+40*f_mod, 0.01, 0.01);

Modulation_filtered = filter(coef,Modulation);
Modulation_filtered_spect = abs(fft(Modulation_filtered)) / length(Modulation_filtered);

figure
subplot(321); plot(t_carrier, Modulation); ylim([-0.5, 1.5]), title('PWM + Carrier')
h4 = subplot(323); plot(f_modu, 2*Modulation_spect(1:length(Modulation_spect)/2+1)); title('PWM + Carrier Spectrum')
h5 = subplot(325); plot(f_modu, 2*Modulation_filtered_spect(1:length(Modulation_filtered_spect)/2+1));...
    title('PWM + Carrier Spectrum Filtered')
linkaxes([h4,h5],'x');

salida = amdemod(Modulation_filtered, fc,f_clk);
subplot(3,2,[2 4 6]); plot(t_carrier, salida-mean(salida)); , title('SALIDA')