

%==================================================================
%This Script produces results in TABLE I in our paper.

%==================================================================

clear
clc

m=20;

n=40:20:120;  %varing n


iter =1000; % number of iterations
rho=0.2; %parameter rho in DCACO

dd1=0.8; dd2=0.95;%parameter dd1 and dd2 in Elad

rho1=0.2;%parameter rho in DMCM-PG


q=10;%the number of initial matrices when fixing n

initiA_mu=zeros(q,1);
time=zeros(q,4);
out_mu=zeros(q,4);


RESULT_times_n=zeros(length(n),4);
RESULT_n=zeros(length(n),6);

for jjj=1:1:length(n)
    
    wb=sqrt((n(jjj)-m)/(m*(n(jjj)-1))); % WB
    
    for j=1:1:q
        A = randn(m,n(jjj));
        A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
        G=A'*A;
        initiA_mu(j) = max(max(abs(G) - eye(n(jjj))));
        
        [time(j,:),out_mu(j,:)]=fourmethods_mu(A,rho,rho1,dd1,dd2,iter);
        
    end
    RESULT_times_n(jjj,1)=mean(time(:,1));
    RESULT_times_n(jjj,2)=mean(time(:,2));
    RESULT_times_n(jjj,3)=mean(time(:,3));
    RESULT_times_n(jjj,4)=mean(time(:,4));
    
    RESULT_n(jjj,1)=wb;
    RESULT_n(jjj,2)=mean(initiA_mu);
    RESULT_n(jjj,3)=mean(out_mu(:,1));
    RESULT_n(jjj,4)=mean(out_mu(:,2));
    RESULT_n(jjj,5)=mean(out_mu(:,3));
    RESULT_n(jjj,6)=mean(out_mu(:,4));
    
end

save('RESULT_times_n.mat','RESULT_times_n')

save('RESULT_n.mat','RESULT_n')
