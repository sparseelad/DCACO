
% ======================================================
% This Script tests DCACO.
% 
% 
% Good luck.
% ======================================================

clear
clc

%%
m = 20; % dimension of the space      
n = 40; % number of vectors


app_wb=1/sqrt(m); % approximate WB
wb=sqrt((n-m)/(m*(n-1))); % WB

iter = 1000; % number of iterations
rho=0.2; %parameter rho

%% initializatio
A = randn(m, n);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));

%% improve initialization by SVD
% [U, ~, V] = svd(A);
% A = U*[eye(m) zeros(m, n-m)]*V';
% A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
% G=A'*A;
% mu_A = max(max(abs(G) - eye(n)))

%% main iteration

[mu_A,mu_newA,newA]=DCACO(A,rho,iter);

%% show results
outs=[app_wb,wb,mu_A,mu_newA];
fprintf('app_wb=%12.8f\nwb=%12.8f\nmu_A=%12.8f\nmu_newA=%12.8f\n',outs);





