clc
clear all
close all

SNR_dB = 10; 
SNR_linear = 10.^(SNR_dB/10.);
N_iter=50000; 
sq2=sqrt(0.5); 
grps = {'b-'; 'r:';'m--';'g-.'};

CDF = zeros(4,50);

for nT = 1:4
    
    nR = nT;
    I = eye(nT);
    
    C = zeros(1,N_iter);
    
    for iter=1:N_iter
        
        H = sq2*(randn(nR,nT)+1j*randn(nR,nT));
        C(iter) = log2(real(det(I+SNR_linear/nT*(H'*H))));
        
    end
    
    [PDF,Rate] = hist(C,50);
    PDF = PDF/N_iter;
    
    for i=1:50
        CDF(nT,i) = sum(PDF(1:i));
    end
    
    plot(Rate,CDF(nT,:),grps{nT},'linewidth',3); 
    hold on
    
end

xlabel('Channel capacity[bps/Hz]'); 
ylabel('CDF')

axis([1 18 0 1]); 
grid on; 
set(gca,'fontsize',10);
legend('{\it N_T}={\it N_R}=1','{\it N_T}={\it N_R}=2',...
    '{\it N_T}={\it N_R}=3','{\it N_T}={\it N_R}=4','location','southeast');
