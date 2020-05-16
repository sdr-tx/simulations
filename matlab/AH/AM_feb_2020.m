clc; ; close all; clear;
%%
BITS = 4;
PWM_STEPS = 2^BITS;
UPSAMPLING = 100;
DATA_CYC = 2;

F_CARRY = 1000e3;                T_CARRY = 1/F_CARRY;

F_CLK = F_CARRY * UPSAMPLING;   T_CLK = 1/F_CLK;

T_PWM = T_CARRY * PWM_STEPS;    F_PWM = 1/T_PWM;

F_MSG = 1e3;                    T_MSG = 1/F_MSG;

N = DATA_CYC * (T_MSG/T_CLK);

%%
t = (0:N-1) * T_CLK;

% RF carrier from PLL
carrier = (square(F_CARRY*2*pi*t) + 1)/2;

msg = 0.5*cos(F_MSG*2*pi*t)+0.5;

DATA_QUANT = 2^BITS/max(msg);

msg_digital = uint8(round(DATA_QUANT * msg));

pwm = uint8(zeros(1,N));

for i = 1:PWM_STEPS:N
    for j = 0:PWM_STEPS-1
        if j < msg_digital(i)
            pwm(i+j) = 1;
        else
            pwm(i+j) = 0;
        end
    end
end

pwm = pwm(1:N);

h1 = subplot(311), plot(t,msg), ylim([-0.2 1.2])
h2 = subplot(312), plot(t,msg_digital), ylim([-PWM_STEPS*0.1 PWM_STEPS*1.1])
h3 = subplot(313), plot(t,pwm), ylim([-0.2 1.2])
linkaxes([h1,h2,h3],'x');

%%

f = (0:N-1)*(F_CLK/N);

carrier_spec = abs(fft(carrier))/N;
msg_spec = abs(fft(msg))/N;
msg_digital_spec = abs(fft(double(msg_digital)))/N;
pwm_spec = abs(fft(double(pwm)))/N;

figure(2);
% subplot(221);   plot(f(1:N/2), mag2db(carrier_spec(1:N/2)));    title('carrier');
% subplot(222);   plot(f(1:N/2), mag2db(msg_spec(1:N/2)));    title('data - analog');
% subplot(223);   plot(f(1:N/2), mag2db(msg_digital_spec(1:N/2)));         title('data - digital');
% subplot(224);   plot(f(1:N/2), mag2db(pwm_spec(1:N/2)));              title('PWM');

subplot(221);   plot(f(1:N/2), carrier_spec(1:N/2));    title('carrier');
subplot(222);   plot(f(1:N/2), msg_spec(1:N/2));    title('data - analog');
subplot(223);   plot(f(1:N/2), (msg_digital_spec(1:N/2)));         title('data - digital');
subplot(224);   plot(f(1:N/2), (pwm_spec(1:N/2)));              title('PWM');
