function [t, y] = cuadrada( Ncyc, duty,  freq_sig, fs)
    T = 1/freq_sig;     %Periodo
    tf = T*Ncyc;        %Periodo de tiempo
    dt = 1/fs;
    t = 0:dt:tf-dt;       %Vector de tiempo
        
    y = zeros(1,fs/freq_sig);
    
    for k=1:T*fs*duty
        y(k) = 1;
    end
    
    
    y = repmat(y, [1 Ncyc]);
    
end

