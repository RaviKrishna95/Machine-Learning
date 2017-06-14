clc
clear all 
close all

SNr_dB = 0:5:20;
SNr_linear = 10.^(SNr_dB/10.);

% Set the number of Iterations.
N_iter = 1000;

N_SNr = length(SNr_dB);
% sq2 = sqrt(0.5);

% Preallocation for channel capacity
C_det = zeros(4,N_SNr);

% Simulating for differeNt number of transmit and receiver aNtennas

for Nt = 2:5
    
    Nr = Nt;
    I = eye(Nr);
    
    for iter=1:N_iter
        
        H = randn(Nr,Nt)+1j*randn(Nr,Nt);
        tmp1 = H'*H/Nt; 
          
        for i=1:N_SNr 
            C_det(Nt-1,i) = C_det(Nt-1,i) + log2(det(I+SNr_linear(i)*tmp1));  
        end
        
    end
    
end

C_det = real(C_det)/N_iter;

% Plotting

plot(SNr_dB,C_det,'-o','linewidth',2)
grid on
xlabel('SNR(dB)','fontsize',10)
ylabel('Channel Capacity (bps/Hz)','fontsize',10)
title('Deterministic MIMO Channel Capacity in Terms of SNR','fontsize',12)
legend('N_T = 2, N_R = 2','N_T = 3, N_R = 3',...
    'N_T = 4, N_R = 4','N_T = 5, N_R = 5','location','Northwest')
