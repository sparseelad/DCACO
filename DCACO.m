function [mu_A,mu_newA,newA]=DCACO(A,rho,iter)

%==================================================================
%This Script implements DCACO described in Algorithm 2 in our paper.

%===================================================================


[m,n]=size(A);
A=bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));

for k=1:1:iter
    Y=sub_diff(G,m);
    W=Y-eye(n);
    w=smatrix_to_vec(W); % can also use Matlab command: w=smatrix_to_vec_simple(W); 
    u=touying_L1ball(w,1/(2*rho));
    z=w-u;
    G=vec_to_smatrix(z,n);
end

newA=incomatrix(G,m);
mu_newA=compu_mu(newA);
end
