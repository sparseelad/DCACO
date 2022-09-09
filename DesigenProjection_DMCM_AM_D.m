function newP=DesigenProjection_DMCM_AM_D(P,D,D_typ,rho,bate,iter)

%======================================================================
% This Script is directly from Algorithm 3 in the paper: 
% "Optimized  Projections for Compressed Sensing via Direct Mutual 
% Coherence Minimization," 
% Signal Process., vol. 151, pp. 45-55, 2018.
 
% 
% Note that: For a fair comparison, we do not implement Algorithm 4 
% in the above paper. In fact, Algorithm 4 in the above paper uses a 
% so-called "continuation trick" to accelerate Algorithm 4 in the 
% above paper
%

%======================================================================

% parameter: rho, bate

%% DMCM_AM_D




A=P*D;

[m,n]=size(A);

A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
% mu_PD=max(max(abs(G) - eye(n)));

a=0.99*rho;
a=1/a;
b=1/bate;
for k=1:1:iter
    
    W=(A'*A-eye(n))*(1/rho);
    v=touying_L1ball(W(:),1);
    V=reshape(v,n,n);
    Z=(1/(a+b))*(a*A+b*P*D-A*(V+V'));
    
    A = bsxfun(@rdivide, Z, sqrt(sum(Z.^2)));
    
    if D_typ==1  
       newP=A;
    else
       newP=A*pinv(D);
    end
    
    
    %newP=A*pinv(D);
    P=newP;
    
    
end
% A=P*D;
% A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
% G=A'*A;
% mu_newPD=max(max(abs(G) - eye(n)));
end