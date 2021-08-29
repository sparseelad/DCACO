
%==================================================================
%This Script produces results in TABLE II in our paper.

%==================================================================




clear
clc

m=20:20:100; %varing m

n=150;


iter =1000; % number of iterations
rho=0.2; %parameter rho in DCACO

dd1=0.8; dd2=0.95;%parameter dd1 and dd2 in Elad

rho1=0.2;%parameter rho in DMCM-PG


q=10; %the number of initial matrices when fixing m

initiA_mu=zeros(q,1);
time=zeros(q,4);
out_mu=zeros(q,4);


RESULT_times_m=zeros(length(m),4);
RESULT_m=zeros(length(m),6);

for jjj=1:1:length(m)
    
    wb=sqrt((n-m(jjj))/(m(jjj)*(n-1))); % WB
    
    for j=1:1:q
        A = randn(m(jjj),n);
        A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
        G=A'*A;
        initiA_mu(j) = max(max(abs(G) - eye(n)));
        
        [time(j,:),out_mu(j,:)]=fourmethods_mu(A,rho,rho1,dd1,dd2,iter);
        
    end
    RESULT_times_m(jjj,1)=mean(time(:,1));
    RESULT_times_m(jjj,2)=mean(time(:,2));
    RESULT_times_m(jjj,3)=mean(time(:,3));
    RESULT_times_m(jjj,4)=mean(time(:,4));
    
    RESULT_m(jjj,1)=wb;
    RESULT_m(jjj,2)=mean(initiA_mu);
    RESULT_m(jjj,3)=mean(out_mu(:,1));
    RESULT_m(jjj,4)=mean(out_mu(:,2));
    RESULT_m(jjj,5)=mean(out_mu(:,3));
    RESULT_m(jjj,6)=mean(out_mu(:,4));
    
end

save('RESULT_times_m.mat','RESULT_times_m')
save('RESULT_m.mat','RESULT_m')