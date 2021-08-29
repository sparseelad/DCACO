

%======================================================
%This Script produces Fig. 3 and Fig. 4 in our paper.

%======================================================

clear
clc

% m = 20; % dimension of the space
% n = 40; % number of vectors
m = 30; % dimension of the space
n = 60; % number of vectors

app_wb = 1/sqrt(m); % performance limits
wb = sqrt((n-m)/(m*(n-1))); %WB

iter = 2000;% number of iterations
rho=0.2; %paramter rho
ss=100; % number of signals

res=zeros(ss,1);

for K=1:1:ss
    
    
    A = randn(m, n);
    A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
    
    % improve initialization by SVD
    
%     [U, ~, V] = svd(A);
%     A = U*[eye(m) zeros(m, n-m)]*V';
%     A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
    
    G=A'*A;
    
    for k=1:1:iter

        Y=sub_diff(G,m);
        W=Y-eye(n);
        w=smatrix_to_vec(W);
        u=touying_L1ball(w,1/(2*rho));
        z=w-u;
        G=vec_to_smatrix(z,n);

    end
    
    
    newA=incomatrix(G,m);
    mu_newA=compu_mu(newA);
    res(K)=mu_newA;
    
end

figure(1)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(1:1:ss,res,'b.','MarkerSize',14);hold on
plot(1:1:ss,mean(res)*ones(ss,1),'m--','LineWidth',3);hold on
plot(1:1:ss,app_wb*ones(ss,1),'m.','MarkerSize',14);hold on
plot(1:1:ss,wb*ones(ss,1),'r.','MarkerSize',14);hold on


xlabel('Number of signals','FontSize',14);
ylabel('Value of \mu','FontSize',14);
h=legend('DCACO','Mean','$1/\sqrt{m}$','Welch bound');
set(h,'interpreter','latex');