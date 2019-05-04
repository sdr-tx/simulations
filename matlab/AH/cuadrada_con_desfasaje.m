function [y, t] = cuadrada_con_desfasaje(Ncyc, Amp, desfasaje, duty, freq_sig, fs)
    %   Ncyc:       Numero de ciclos que se quiere de la señal
    %   Amp:        Amplitud de la señal
    %   desfasaje:  Desfasaje de la señal: desfasaje 0 comienza en 1 T/2
    %   duty:       Duty cicle (50% = 0.5)
    %   freq_sig:   Frecuencia de la señal que se quiere sacar
    %   fs:         frecuencia de muestreo de la señal (fs = 1/dt);
    
    T = 1/freq_sig;     %Periodo
    tf = T*Ncyc;        %Periodo de tiempo
    dt = 1/fs;
    t = 0:dt:tf-dt;       %Vector de tiempo
    phase_n = -desfasaje * fs/freq_sig;
        
    y = zeros(1,fs/freq_sig);
    
    for k=1:T*fs*duty
        y(k) = 1;
    end
    
    y = circshift(y',phase_n);
    y = Amp .* y';
    y = repmat(y, [1 Ncyc]);
 
    
end

