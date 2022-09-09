

%==================================================================
%This Script produces results in TABLE 1 in our paper.

%==================================================================

clear
clc

m=20;

n=[40,100,500];  %varing n


iter =1000; % number of iterations
rho=0.2; %parameter rho in DCACO

dd1=0.8; dd2=0.95;%parameter dd1 and dd2 in Elad

alpha=0.5;%parameter alpha in ISPM

rho1=0.2;%parameter rho in DMCM-PG


q=5;%the number of initial matrices when fixing n

initiA_mu=zeros(q,1);
time=zeros(q,5);
out_mu=zeros(q,5);


RESULT_times_n=zeros(length(n),5);
RESULT_n=zeros(length(n),7);

for jjj=1:1:length(n)
    
    wb=sqrt((n(jjj)-m)/(m*(n(jjj)-1))); % WB
    
%     a_mn=3*n(jjj)-m^2-2*m;
%     b_mn=(m+2)*(n(jjj)-m);
%     lb=sqrt(a_mn/b_mn); % LB
%     desired_coh=1.1*max(wb,lb); %parameter desired_coh in ISPM
    desired_coh=1.1*wb;
    for j=1:1:q
        A = randn(m,n(jjj));
        A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
        G=A'*A;
        initiA_mu(j) = max(max(abs(G) - eye(n(jjj))));
        
        [time(j,:),out_mu(j,:)]=fivemethods_mu(A,rho,rho1,dd1,dd2,desired_coh,alpha,iter);
        
    end
    RESULT_times_n(jjj,1)=mean(time(:,1));
    RESULT_times_n(jjj,2)=mean(time(:,2));
    RESULT_times_n(jjj,3)=mean(time(:,3));
    RESULT_times_n(jjj,4)=mean(time(:,4));
    RESULT_times_n(jjj,5)=mean(time(:,5));
    
    RESULT_n(jjj,1)=wb;
    RESULT_n(jjj,2)=mean(initiA_mu);
    RESULT_n(jjj,3)=mean(out_mu(:,1));
    RESULT_n(jjj,4)=mean(out_mu(:,2));
    RESULT_n(jjj,5)=mean(out_mu(:,3));
    RESULT_n(jjj,6)=mean(out_mu(:,4));
    RESULT_n(jjj,7)=mean(out_mu(:,5));
    
end

save('RESULT_times_n.mat','RESULT_times_n')

save('RESULT_n.mat','RESULT_n')
