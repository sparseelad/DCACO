function []=optimal_parameter_Fig()

%====================================================
%This Script produces Fig. 1 and Fig. 2 in our paper.

%====================================================

% m = 20;% dimension of the space
% n = 40;% number of vectors
m = 30;% dimension of the space
n = 60;% number of vectors
% initializatio
A = randn(m, n);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));

% improve initialization by SVD
% [U, ~, V] = svd(A);
% A = U*[eye(m) zeros(m, n-m)]*V';
% A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
% G=A'*A;
% mcA = max(max(abs(G) - eye(n)))

app_wb=1/sqrt(m); % performance limits
wb=sqrt((n-m)/(m*(n-1))); %WB

iter = 2000; % number of iterations

res=zeros(iter+1,13);
count=1;
for rho=[0.2,0.4,0.6,0.8,1,2,4,6,8,10,20,50,100]
    
    [res(:,count)]=DCACO(A,rho,iter);
    count=count+1;
end

figure(1)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(0:1:iter,res(:,1),'b-','LineWidth',3);hold on
plot(0:1:iter,res(:,2),'r-','LineWidth',3);hold on
plot(0:1:iter,res(:,3),'g-','LineWidth',3);hold on
plot(0:1:iter,res(:,4),'c-','LineWidth',3);hold on
plot(0:1:iter,res(:,5),'k-','LineWidth',3);hold on
plot(0:1:iter,wb*ones(iter+1,1),'r--','LineWidth',3);hold on
axis([0 2000 0 1]);
xlabel('Number of iterations','FontSize',14);
ylabel('Value of \mu','FontSize',14);
legend('\rho=0.2','\rho=0.4','\rho=0.6','\rho=0.8','\rho=1','Welch bound')

figure(2)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(0:1:iter,res(:,5),'b-','LineWidth',3);hold on
plot(0:1:iter,res(:,6),'r-','LineWidth',3);hold on
plot(0:1:iter,res(:,7),'g-','LineWidth',3);hold on
plot(0:1:iter,res(:,8),'c-','LineWidth',3);hold on
plot(0:1:iter,res(:,9),'k-','LineWidth',3);hold on
plot(0:1:iter,wb*ones(iter+1,1),'r--','LineWidth',3);hold on
axis([0 2000 0 1]);
xlabel('Number of iterations','FontSize',14);
ylabel('Value of \mu','FontSize',14);
legend('\rho=1','\rho=2','\rho=4','\rho=6','\rho=8','Welch bound')


figure(3)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(0:1:iter,res(:,9),'b-','LineWidth',3);hold on
plot(0:1:iter,res(:,10),'r-','LineWidth',3);hold on
plot(0:1:iter,res(:,11),'g-','LineWidth',3);hold on
plot(0:1:iter,res(:,12),'c-','LineWidth',3);hold on
plot(0:1:iter,res(:,13),'k-','LineWidth',3);hold on
plot(0:1:iter,wb*ones(iter+1,1),'r--','LineWidth',3);hold on
axis([0 2000 0 1]);
xlabel('Number of iterations','FontSize',14);
ylabel('Value of \mu','FontSize',14);
legend('\rho=8','\rho=10','\rho=20','\rho=50','\rho=100','Welch bound')



end



function [out]=DCACO(A,rho,iter)

[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
out=zeros(iter+1,1);
out(1)=mu_A;

for k=1:1:iter
    Y=sub_diff(G,m);
    W=Y-eye(n);
    w=smatrix_to_vec(W);
    u=touying_L1ball(w,1/(2*rho));
    z=w-u;
    G=vec_to_smatrix(z,n);
    
    
    newA=incomatrix(G,m);
    mu_newA=compu_mu(newA);
    
    out(k+1)=mu_newA;
    
end

end
