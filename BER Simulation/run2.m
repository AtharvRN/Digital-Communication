
function run2()
    clc;
    clear all;
    M_arr = [2 4 16 64 256];
    %M_arr=[4];
    Nf = 1000;
    Na = 1024;
   for ii = 1:length(M_arr)
     M = M_arr(ii);
    disp(sprintf('M = %d',M));
    figure(ii);
    clf(figure(ii));
    PAM_sim(M,ii,Nf,Na);
    PSK_sim(M,ii,Nf,Na);
    QAM_sim(M,ii,Nf,Na);
   end
end

function PAM_sim(M,ii,Nf,Na)
    disp('PAM');
    P = zeros(1,9);
    P_theoretical = zeros(1,9);
    % No of bits
    
    m = log2(M);
   Eb = 1;
   Es = m*Eb;
   % Generate a LUT
     D =  sqrt(3*Es/(M^2-1));
    lookup = D.*[-(M-1):2:(M-1)];
    % Gray table generation
     gray = zeros(1,M);
    for tt = 1:M
        gray(tt) = bitxor(tt-1,floor((tt-1)/2));
    end
    
    SNR = 0;
    for count = 1:1:9
    
        Ne_total = 0;
        %Nf = 100;
        % Unit Bit Energy 
        
        N0 = Es/power(10,SNR/10);
        std_dev = sqrt(N0/2);
        mean = 0;
        for times = 1:Nf
            %Na = 1000;
            % a = 0 or 1
            a = rand(1,Na)>0.5 ;
            
            v = [];
            w = [];
            r= [];
           
            for i=1:Na/m
                
                x = a((i-1)*m +1:i*m);
                %disp(x);
                dec= get_decimal(x);

                % PAM MAPPING
                v(i) = D*lookup(gray(dec+1)+1);
                %s disp(v(i));
             
                % Noise
                w(i) = std_dev*randn + mean;
                %disp(w(i));
            % Transmitted Signal
                r(i) = v(i) + w(i);
                %disp(r(i));

            
        [dem, index] = min(abs(lookup-r(i)));
        
         recv = find(gray == index-1);
         detected = get_bits(recv-1,m);
        
         for j =1:m
              if(detected(j) ~= x(j))
                   Ne_total = Ne_total+1;
              end
         end

            end
        end

        
        P(count) = Ne_total/(Na*Nf);
        P_theoretical(count) = (2*((M-1)/M)*qfunc(sqrt((6*m*Es)/((M^2-1)*N0))))/m;
        X = sprintf('SNR =  %d dB Practical BER PAM = %d \n',SNR,P(count));

        Y = sprintf('SNR = %d dB  Theoretical BER  PAM= %d \n',SNR,P_theoretical(count));
         
        
        SNR = SNR+ 1;
        format long
        disp(X);
        disp(Y);
    end
    
    
    clf(figure(ii));
    x = [0:8];
    semilogy(x,P);
    hold on;
    semilogy(x,P_theoretical);
    hold on;
    
    title(sprintf('BER vs SNR for M = %d',M));
    xlabel('SNR (dB)');
    ylabel('BER');
    legend('Simulation BER PAM','Theoretical BER PAM');
    
    
end


function PSK_sim(M,ii,Nf,Na)
    P = zeros(1,9);
    P_theoretical = zeros(1,9);
    disp('PSK');
    SNR = 0;
    
  
    % No of bits
     m = log2(M);
     Eb = 1;
     Es = m*Eb;
     D =  sqrt(Es); 
                lookup = D*exp((1i*pi/M)*([1:2:2*M-1]));
                % Gray table generation
                gray = [];
                for tt = 1:M
                    gray(tt) = bitxor(tt-1,floor((tt-1)/2));
                end
                gray;
    for count = 1:1:9
    
        Ne_total = 0;
        %Nf = 100;
        for times = 1:Nf
            %Na = 1000;
            % a = 0 or 1
            a = rand(1,Na)>0.5 ;
            % Bit Energy 
            v = [];
            w = [];
            r= [];
            % M-PAM
            
            N0 = Es/power(10,SNR/10);
            std_dev = sqrt(N0/2);
            mean = 0;
            % Generate a LUT
                
               
            for i=1:Na/m
                
                x = a((i-1)*m +1:i*m);
                %disp(x);
                dec= get_decimal(x);

                % PSK MAPPING
                v(i) = D*lookup(gray(dec+1)+1);
                %disp(v(i));
             
                % Noise
                w(i) = sqrt(1/2)*(std_dev*randn + 1i*std_dev*randn);
                %disp(w(i));
            % Transmitted Signal
                r(i) = v(i) + w(i);
                %disp(r(i));

            
        [dem, index] = min(abs(lookup-r(i)));
        
         recv = find(gray == index-1);
         detected = get_bits(recv-1,m);
        
                for j =1:m
                    if(detected(j) ~= x(j))
                        Ne_total = Ne_total+1;
                    end
                end

            end
        end

        
        P(count) = Ne_total/(Na*Nf);
        P_theoretical(count) = ((M-1)/m)*qfunc(D*sin(pi/M)*sqrt(2/N0));
        %P_theoretical(count) = qfunc(sqrt(2/N0));
        X = sprintf('SNR =  %d dB Practical BER PSK = %d \n',SNR,P(count));

        Y = sprintf('SNR = %d dB  Theoretical BER PSK = %d \n',SNR,P_theoretical(count));
         
        
        SNR = SNR+ 1;
        format long
        disp(X);
        disp(Y);
    end
    
    %figure(ii);
  
    x = [0:8];
    semilogy(x,P);
    hold on;
    semilogy(x,P_theoretical);
    hold on;
    
    title(sprintf('BER vs SNR for M = %d',M));
    xlabel('SNR (dB)');
    ylabel('BER');
    legend('Simulation BER PAM','Theoretical BER PAM','Simulation BER  PSK','Theoretical BER PSK');
    

end

function QAM_sim(M,ii,Nf,Na)

    P = zeros(1,9);
    P_theoretical = zeros(1,9);
    disp('QAM');
    SNR = 0;
    
  
    % No of bits
     m = log2(M);
     SNR = 0;
     Eb  =1 ;
     Es = m*Eb;
     lookup = zeros(1,M);
     gray = zeros(1,M);
     
      D =  sqrt(1.5*Es/(M-1));
                k = sqrt(M);
       for n = 0:M-1
          lookup(n+1) = D*((2*mod(n,k) + 1-k)*((-1)^(floor(n/k)))+1i*(2*floor(n/k)+1-k));
       end
                % Gray table generation
               
        for tt = 1:M
              gray(tt) = bitxor(tt-1,floor((tt-1)/2));
        end
       
    for count = 1:1:9
    
        Ne_total = 0;
        %Nf = 100;
        for times = 1:Nf
            %Na = 1000;
            % a = 0 or 1
            a = rand(1,Na)>0.5 ;
            % Bit Energy 
            Eb = 1;
            v = [];
            w = [];
            r= [];
            % M-PAM
            
            N0 = Eb/power(10,SNR/10);
            std_dev = sqrt(N0/2);
            mean = 0;
            % Generate a LUT
                
                
            
               
            for i=1:Na/m
                
                x = a((i-1)*m +1:i*m);
                %disp(x);
                dec= get_decimal(x);

                % PSK MAPPING
                v(i) = D*lookup(gray(dec+1)+1);
               %s disp(v(i));
             
                % Noise
                w(i) = (std_dev*randn + 1i*std_dev*randn);
                %disp(w(i));
            % Transmitted Signal
                r(i) = v(i) + w(i);
                %disp(r(i));

            
        [dem, index] = min(abs(lookup-r(i)));
        
         recv = find(gray == index-1);
         detected = get_bits(recv-1,m);
        
                for j =1:m
                    if(detected(j) ~= x(j))
                        Ne_total = Ne_total+1;
                    end
                end

            end
        end

        
        P(count) = Ne_total/(Na*Nf);
        P_theoretical(count) = 4*(k-1)/(m*k)*qfunc((1.5/(M-1))*log2(M)*sqrt(2/N0));
        X = sprintf('SNR =  %d dB Practical QAM  = %d \n',SNR,P(count));

        Y = sprintf('SNR = %d dB  Theoretical BER QAM = %d \n',SNR,P_theoretical(count));
         
        
        SNR = SNR+ 1;
        format long
        disp(X);
        disp(Y);
    end
    
    %figure(ii);
  
    x = [0:8];
    semilogy(x,P);
    hold on;
    semilogy(x,P_theoretical);
    hold on;
    
    %title(sprintf('BER vs SNR for M = %d',M));
    %xlabel('SNR (dB)');
    %ylabel('BER');
    legend('Simulation BER PAM','Theoretical BER PAM','Simulation BER  PSK','Theoretical BER PSK','Simulation BER  QAM','Theoretical BER QAM');
    

end





    function level = get_decimal(x)
        m = length(x);
        level = 0;
    for i = m:-1:1
        if x(i) == 1 
            level = level + 2^(m-i);
        end
    end
    end

    function bit = get_bits(num,m)
        bit = zeros(1,m);
        temp = num;
        for i = m:-1:1
        if(mod(temp,2) == 1)
            bit(i) = 1;
        end
        temp = fix(temp/2);
        end
    end
        


   



