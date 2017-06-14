clc
clear all 
close all

Rt = [1 0.76*exp(0.17j*pi) 0.43*exp(0.35j*pi) 0.25*exp(0.53j*pi); ...
0.76*exp(-0.17j*pi) 1 0.76*exp(0.17j*pi) 0.43*exp(0.35j*pi); ...
0.43*exp(-0.35j*pi) 0.76*exp(-0.17j*pi) 1 0.76*exp(0.17j*pi); ...
0.25*exp(-0.53j*pi) 0.43*exp(-0.35j*pi) 0.76*exp(-0.17j*pi) 1];

Rr = [1 0.76*exp(0.17j*pi) 0.43*exp(0.35j*pi) 0.25*exp(0.53j*pi); ...
0.76*exp(0.17j*pi) 1 0.76*exp(0.17j*pi) 0.43*exp(0.35j*pi); ...
0.43*exp(0.35j*pi) 0.76*exp(0.17j*pi) 1 0.76*exp(0.17j*pi); ...
0.25*exp(0.53j*pi) 0.43*exp(0.35j*pi) 0.76*exp(0.17j*pi) 1];

% Vary the Input SNr

SNr_dB = 0:5:27;
SNr_linear = 10.^(SNr_dB/10.);

% Set the number of Iterations.
N_iter = 1000;

N_SNr = length(SNr_dB);
sq2 = sqrt(0.5);

% Preallocation for channel capacity
C_iid = zeros(2,N_SNr);
C_corr = zeros(2,N_SNr);
    
% Simulating for differeNt number of transmit and receiver aNtennas

for Nt = 3:4
    
    Nr = Nt;
    I = eye(Nr);
    
    for iter=1:N_iter
        
        H_iid = sq2*(randn(Nr,Nt)+1j*randn(Nr,Nt));
        H_corr = (Rr(1:Nt,1:Nt)^(1/2))*H_iid*(Rt(1:Nt,1:Nt)^(1/2));
        tmp1 = H_iid'*H_iid/Nt; 
        tmp2 = H_corr'*H_corr/Nt;
        
        for i=1:N_SNr 
            C_iid(Nt-2,i) = C_iid(Nt-2,i) + log2(det(I+SNr_linear(i)*tmp1));  
            C_corr(Nt-2,i) = C_corr(Nt-2,i) + log2(det(I+SNr_linear(i)*tmp2));
        end
        
    end
    
end

C_iid = real(C_iid)/N_iter;
C_corr = real(C_corr)/N_iter;

% Plotting

plot(SNr_dB,[C_iid;C_corr],'-o','linewidth',2)
grid on
xlabel('SNR(dB)','fontsize',12)
ylabel('Channel Capacity (bps/Hz)','fontsize',12)
title('Capacity of iid and correlated channels ','fontsize',14)
legend('3X3 iid channel','4X4 iid channel',...
    '3X3 correlated channel','4X4 correlated channel','location','Northwest')



