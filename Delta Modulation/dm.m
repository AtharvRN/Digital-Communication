clc;
clear all;
%Sampling Frequency
fs = 100;
Ts = 1/fs;

T = 2;
numSamples = T/Ts; 
%defining a sine wave with appropriate paramters
t = linspace(0,T,numSamples);
f = 1;
m = sin(2*pi*f*t);
figure(1)
plot(m);
title('Sinusoidal Signal');
hold on
%% TRADIONAL DELTA MODULATION
%quantized signal
xq =[0];
y1 =[0];
%step size
delta = 0.05;
for i = 1:numSamples-1
    if m(i) >= xq(i)
        xq(i+1) = xq(i) + delta;
        y1(i) = 1;
    else 
        xq(i+1) = xq(i) - delta;
        y1(i) = 0;
    end
end
% To obtain relevant results, uncomment the below lines
%stairs(xq);
%hold off;
%plot(y1);

%% ADAPTIVE DELTA MODULATION
% quantized signal
xqa =[0];
err_prev = 0;
y2 = [0];

% initial step size
delta_0 = 0.1;
delta_prev = delta_0;
for i = 1:numSamples-1
    err = m(i) - xqa(i);

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
title('Transmitted Signal')
legend({'Sinusoidal Signal','Adaptive Delta Modulation'},'Location','southwest')
figure(2);
plot(y2);
title('Bit Stream');
hold off;


%% RECONSTRUCTION FROM BIT STREAM Y2
xr =[0];
err_prev = 0;
% initial step size
delta_i = 0.1;
delta_prev = delta_i;
for i = 1:numSamples-1
  
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

figure(3);
stairs(xr);
%% 
hold on;
w = 0.1*pi;
f = w/(2*pi);
m_ = lowpass(xr,f);

plot(m_);
title('Reciever Side');
legend({'Reconstructed signal from bit stream','After passing through LPF'},'Location','southwest')
