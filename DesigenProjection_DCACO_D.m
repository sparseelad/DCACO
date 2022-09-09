function newP=DesigenProjection_DCACO_D(P,D,D_typ,rho,iter)

% ===========================================================================
% This Script produces an optimized projection matrix by applying DCACO when
% given a dictionary D.

%
% First, applying DCACO with an
% initial effective dictionary A = PD to produce an incoherent
% effective dictionary A_DCACO.
% Second, solving a least squares
% problem, we obtain an optimized projection matrix P_DCACO
% satisfying A_DCACO is approximately equal to P_DCACOD.
% 
% ===========================================================================

%%DCACO
A=P*D;

[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
%mu_PD=max(max(abs(G) - eye(n)));

for k=1:1:iter
    Y=sub_diff(G,m); 
    W=Y-eye(n);
    w=smatrix_to_vec(W);
    u=touying_L1ball(w,1/(2*rho));
    z=w-u;
    G=vec_to_smatrix(z,n);  
   
    
end

newA=incomatrix(G,m);

%%solving a least squares problem
if D_typ==1   
   newP=newA;
else
   newP=newA*pinv(D);
end

%newP=newA*pinv(D);


end
