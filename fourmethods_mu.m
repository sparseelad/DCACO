function [time,out_mu]=fourmethods_mu(A,rho,rho1,dd1,dd2,iter)

% =================================================================
% This Script computes the mutual coherence of a matrix A 
% by running four solvers, that is, DCACO, Elad, SIDCO and DMCM_PG.
% 
% ==================================================================

time=zeros(1,4);

out_mu=zeros(1,4);


[time(1),out_mu(1),~,~]=DCACO(A,rho,iter);

[time(2),out_mu(2),~,~]=Elad(A,dd1,dd2,iter);

[time(3),out_mu(3),~,~]=SIDCO(A, iter);
[time(4),out_mu(4),~,~]=DMCM_PG(A,rho1,iter);


end




function [time,mu_newA,k,newG]=DCACO(A,rho,iter)

% =====================================================
% also see the Script of DCACO in our software package.
% ======================================================




tic;
[m,n]=size(A);
A=bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));

% out=zeros(iter+1,1);
% out(1)=mu_A;

for k=1:1:iter
    Y=sub_diff(G,m);
    W=Y-eye(n);
    w=smatrix_to_vec(W);
    u=touying_L1ball(w,1/(2*rho));
    z=w-u;
    G=vec_to_smatrix(z,n);
    
    newA=incomatrix(G,m);
    mu_newA=compu_mu(newA);
    %out(k+1)=mu_newA;
    
    %stop condition
    if (abs(mu_newA-mu_A)/abs(mu_A)<1e-6)
        break;
    end
    mu_A=mu_newA;
    
    
    
end
newA=incomatrix(G,m);
newA = bsxfun(@rdivide, newA, sqrt(sum(newA.^2)));
newG=newA'*newA;

time = toc;

%fprintf('time_DCACO=%12.8f\n mu_DCACO=%12.8f\n iter_DCACO=%6i\n',[time,mu_newA,k]);

end


function [time,mu_newA,k,newG]=Elad(A,dd1,dd2,iter)

%======================================================================
% This Script is directly from Fig. 2.2 in Chapter 2 in 
% M. Elad's book: "Sparse and Redundant Representations: From Theory to
% Applications in Signal and Image Processing," Springer, New York, 2010.
% 
% This Script is also indirectly from the software package: "SparseLab
% \SparseLab2.1-Core\CompSense" which can be downloaded free from the 
% website: https://elad.cs.technion.ac.il/software/
% 
% More details on the main idea of this Script, please see the paper: 
% "Optimized Projections for Compressed Sensing," 
% IEEE Trans. Signal Process., vol. 55, no. 12, pp. 5695-5702, 2007.

%======================================================================

% % parameter
% dd1: Relative number of the off-diagonal Gram matrix entries to shrink.
% dd2: Shrink factor to use.
% 
% 
% % the algorithm of Elad





tic;
methodSVD=2; % SVD is done by the power-method and not directly

[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
% out=zeros(iter+1,1);
% out(1)=mu_A;

%Threshold=dd1;

for k=1:1:iter
    
    
    % shrink the high inner products
    %     pos=find(abs(G(:))>=Threshold & abs(G(:)-1)>1e-6);
    %     G(pos)=G(pos)*dd2;
    
    gg=sort(abs(G(:)));
    pos=find(abs(G(:))>gg(round(dd1*(n^2-n))) & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    
    
    % reduce the rank back to m
    if methodSVD==1,
        [U,S,V]=svds(G,m);
        S=S(1:m,1:m);
        U=U(:,1:m);
    elseif methodSVD==2,
        U=randn(n,m);
        for jjj=1:1:10, U=G*U; U=orth(U); end;
        S=diag(diag(U'*G*U));
    end
    G=U*S*U';
    newA=S.^0.5*U';
    
    
    
    newA = bsxfun(@rdivide, newA, sqrt(sum(newA.^2)));% Normalize the columns
    G=newA'*newA;
    mu_newA=max(max(abs(G) - eye(n)));
    
    %out(k+1)=mu_newA;
    
    %stop condition
    if (abs(mu_newA-mu_A)/abs(mu_A)<1e-6)
        break;
    end
    mu_A=mu_newA;
    
    
end
newG=G;
time = toc;

%fprintf('mu_Elad=%12.8f\n',out(end));
end




function [time,mu,k,newG] = SIDCO(A, iter)

%======================================================================
%This Script is directly from the website:
%https://github.com/cristian-rusu-research/SIDCO

%
% Note that: In the original SIDCO, the well-known CVX, a package for
% specifying and solving convex programs, was used
% to solve a series of convex subproblems in each iteration.
% Because we are unable to run CVX on our computing platform
% and in turn we find that a series of convex subproblems in
% each iteration of SIDCO fall into the scheme of the well-known
% Chambolle-Pock primal-dual method (CP method
% for short), we use the computationally simple and efficient CP
% method to solve these convex subproblems in SIDCO.
%

% More details on CVX, please see the the website: http://cvxr.com/cvx/
% 
% More details on CP method, please see the paper: 
% "A first-order primal-dual algorithm for
% convex problems with applications to imaging," 
% J. Math. Imag. Vis., vol. 40, no. 1, pp. 120-145, 2011.
% 
% The function "PC_sub" in SIDCO is the concrete implementation of the CP
% method to solve these convex subproblems in SIDCO.
% More details on the function "PC_sub" in SIDCO, please see the Script of
% the function "PC_sub" in our software package.
% 
% More details on the main idea of this Script, please see the paper: 
% C. Rusu and N. Gonz\'{a}lez-Prelcic, "Designing Incoherent Frames
% Through Convex Techniques for Optimized Compressed Sensing,"
% IEEE Trans. Signal Process., vol. 64, no. 9, pp. 2334-2344, 2016.
%======================================================================



%% SIDCO




tic;
[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));

out=zeros(iter+1,1);
out(1)=mu_A;

for k=1:1:iter
    ordering = randperm(n);
    
    for i = ordering
        target = A(:, i);
        Awork = A(:, [1:i-1 i+1:end]);
        
        %% change signs to simplify optimization problem
        signs = (Awork'*target <= 0);
        for j=1:n-1
            if (signs(j) == 1)
                Awork(:, j) = -Awork(:,j);
            end
        end
        
        corre = Awork'*target;
        T = max(corre);
        
        
        newAi=PC_sub(Awork',target,sqrt(1-T^2));
        
        newAi = newAi/norm(newAi);
        A(:, i) = newAi;
        
    end
    %% new frame, new mutual coherence
    mu = max(max(abs(A'*A) - eye(n)));
    out(k+1)=mu;
    %%stop condition
    if (abs(out(k+1)-out(k))/abs(out(k))<1e-6)
        break;
    end
    % check convergence
    if (k >= 5)
        if (std(out(end-2:end)) < 10e-6)
            [U, ~, V] = svd(A);
            A = U*[eye(m) zeros(m, n-m)]*V';
            A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
            ordering = randperm(n);
        end
    end
    
end
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
newG=A'*A;

time = toc;

%fprintf('time_SIDCO=%12.8f\n mu_SIDCO=%12.8f\n iter_SIDCO=%6i\n',[time,mu,k]);

end




function [time,mu_newA,k,newG]=DMCM_PG(A,rho1,iter)

%======================================================================
% This Script is directly from Algorithm 1 in the paper: 
% "Optimized  Projections for Compressed Sensing via Direct Mutual 
% Coherence Minimization," 
% Signal Process., vol. 151, pp. 45-55, 2018.
 
% 
% Note that: For a fair comparison, we do not implement Algorithm 2 
% in the above paper. In fact, Algorithm 2 in the above paper uses a 
% so-called "continuation trick" to accelerate Algorithm 1 in the 
% above paper
%

%======================================================================

% parameter: rho

%% DMCM_PG




tic;
[~,n]=size(A);

A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
% out=zeros(iter+1,1);
% out(1)=mu_A;

a=0.99*rho1;
for k=1:1:iter
    
    W=(A'*A-eye(n))*(1/rho1);
    v=touying_L1ball(W(:),1);
    V=reshape(v,n,n);
    Z=A-a*A*(V+V');
    A = bsxfun(@rdivide, Z, sqrt(sum(Z.^2)));
    
    
    mu_newA=max(max(abs(A'*A) - eye(n)));
    
    %out(k+1)=mu_newA;
    %stop condition
    if (abs(mu_newA-mu_A)/abs(mu_A)<1e-6)
        break;
    end
    mu_A=mu_newA;
    
end
newG=A'*A;

time = toc;
%fprintf('mu_DMCM-PG=%12.8f\n',out(end));
end