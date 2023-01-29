clc;
clear all;


%SNR = 10*log10(Eb/N0) dB;


P =[];
P_theoretical = [];
disp('------Transmitting 100 Blocks--------');
for SNR = 0:1:10
    num_block_errors = 0;
    Ne_total = 0;
    Nf = 100;
    for times = 1:Nf
        Na = 1000;
        % a = 0 or 1
        a = rand(1,Na)>0.5 ;
        % Bit Energy 
        Eb = 1;
        v = (1-2*a)*sqrt(Eb);
        w = [];
        r= [];

        N0 = Eb/power(10,SNR/10);
        std_dev = sqrt(N0/2);
        mean = 0;
        %w =  std_dev.*randn(Na,1) + mean;
        a_detected = [];

        block_error  =0;
        for i =1:Na
            % Noise
            w(i) = normrnd(mean,std_dev);
            % Transmitted Signal
            r(i) = v(i) + w(i);

            % Detection
            if r(i) >=0
                a_detected(i) = 0;
             else 
                a_detected(i) = 1;
            end
            % Counting Errors
            if(a_detected(i) ~= a(i))
               
                Ne_total = Ne_total +1;
            end


        end
    end
  
   P(SNR+1) = Ne_total/(Na*Nf);
   X = sprintf('SNR =  %d dB BER  = %d ',SNR,P(SNR+1));

   % Theoretical BER
   P_theoretical(SNR+1) = 0.5*erfc(sqrt(10^(SNR/10)));
   format long
   disp(X);


end
disp('-------------  50 Block Error BER----------');
P_sim =[];
blocks_transmitted = [];
for SNR = 0:1:10
    num_block_errors = 0;
    Ne_total = 0;
    
    block_no = 0;
    % Loop will keep on running till there are atleast 50 block errors
    while true

        Na = 1000;
        a = rand(1,Na)>0.5 ;
        Eb = 1;
        v = (1-2*a)*sqrt(Eb);
        w = [];
        r= [];
        N0 = Eb/power(10,SNR/10);
        std_dev = sqrt(N0/2);
        mean = 0;
        
        a_detected = [];

        block_error  =0;
        for i =1:Na
            w(i) = normrnd(mean,std_dev);
            r(i) = v(i) + w(i);

            if r(i) >=0
                a_detected(i) = 0;
             else 
                a_detected(i) = 1;
            end

            if(a_detected(i) ~= a(i))
               
                Ne_total = Ne_total +1;
                block_error = 1;
            end


        end
        % Checking block error
        if(block_error == 1)
            num_block_errors = num_block_errors +1;
        end
        if(num_block_errors > 50)
            %disp('Number of block error crosses 50');
            break;
        end

        block_no = block_no+1;
    end
  
   P_sim(SNR+1) = Ne_total/(Na*block_no);
   blocks_transmitted(SNR+1) = block_no;
  
   format long
   X = sprintf('Number of transmitted blocks  is %d',block_no);
   disp(X);
   X = sprintf('SNR =  %d dB BER  = %d ',SNR,P_sim(SNR+1));
   disp(X);


end

 x  = 0:1:10;
close all;
figure(1);
semilogy(x,P);
hold on;
semilogy(x,P_theoretical);
hold on;
semilogy(x,P_sim);

title('BER vs SNR');
xlabel('SNR (dB)');
ylabel('BER');
legend('Simulation BER with 100 blocks transmitted','Theoretical BER','BER - 50 Block Errors');

figure(2)
semilogy(x,blocks_transmitted);
title('Number of blocks transmitted vs SNR (for 50 block errors)')
xlabel('SNR (dB)')
ylabel('Blocks transmitted')





