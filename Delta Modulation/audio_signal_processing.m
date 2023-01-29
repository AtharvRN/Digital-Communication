
clc;
clear all;
[y,fs] = audioread('audio_file.wav');
signal =  100*y(10000:10100,1);
%sound(signal,fs);
figure(1);
plot(signal);
title('100 samples of the audio signal sampled at 48000 Hz - Transmitter End')
%fvtool(signal,1);
hold on;

%% ADAPTIVE DELTA MODULATION
% quantized signal
xqa =[0];
err_prev = 0;
y2 = [0];

% initial step size
delta_0 = 0.1;
delta_prev = delta_0;
for i = 1:length(signal)-1
    err = signal(i) - xqa(i);

    if err >=0
        err_curr = 1;
        y2(i) = 1;
    else 
        err_curr = -1;
        y2(i) = 0;
     end;
    delta_ = abs(delta_prev)*err_curr + delta_0*err_prev
    xqa(i+1) = xqa(i) + delta_;
    delta_prev = delta_;
    err_prev = err_curr;
    
end
stairs(xqa);
legend({'Audio signal','After ADM'},'Location','southwest')



%% Reconstruction from y2 Same procedure as earlier
xr =[0];
err_prev = 0;
% initial step size
delta_i = 0.1;
delta_prev = delta_i;
for i = 1:length(signal)-1
  
    if y2(i) == 1
        err_curr = 1;
    else 
        err_curr = -1;
    end
    delta_ = abs(delta_prev)*err_curr + delta_i*err_prev;
    xr(i+1) = xr(i) + delta_;
    delta_prev = delta_;
    err_prev = err_curr;
    
end

figure(2);
stairs(xr);
hold on;
w = 0.1*pi;
f = w/2*pi;
m_ = lowpass(xr,0.1*pi);
plot(m_);
title('Reciever End')
legend({'Reconstructed signal from bit stream','After passing through LPF'},'Location','southwest')




