function [time,out_mu]=fivemethods_mu(A,rho,rho1,dd1,dd2,desired_coh,alpha,iter)

% =================================================================
% This Script computes the mutual coherence of a matrix A 
% by running five solvers, that is, DCACO, Elad, SIDCO£¬ISPM and DMCM_PG.
% 
% ==================================================================

time=zeros(1,5);

out_mu=zeros(1,5);


[time(1),out_mu(1)]=DCACO(A,rho,iter);

[time(2),out_mu(2)]=Elad(A,dd1,dd2,iter);

[time(3),out_mu(3)]=SIDCO(A, iter);

[time(4),out_mu(4)]=ISPM(A, desired_coh, alpha,iter);

[time(5),out_mu(5)]=DMCM_PG(A,rho1,iter);


end




function [time,best_mu]=DCACO(A,rho,iter)

% =====================================================
% also see the Script of DCACO in our software package.
% ======================================================




tic;
[m,n]=size(A);
A=bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));

out=zeros(iter+1,1);
out(1)=mu_A;
err_min=1;

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
    
    
    
    newA = bsxfun(@rdivide, newA, sqrt(sum(newA.^2)));
    newG=newA'*newA;
    mu_newA=max(max(abs(newG) - eye(n)));
    
    out(k+1)=mu_newA;
    
    %% update best
    if out(k+1) < err_min
       err_min = out(k+1);
       best_A = newA;
    end
    
    
end

best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));


time = toc;

%fprintf('time_DCACO=%12.8f\n mu_DCACO=%12.8f\n iter_DCACO=%6i\n',[time,mu_newA,k]);

end


function [time,best_mu]=Elad(A,dd1,dd2,iter)

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
out=zeros(iter+1,1);
out(1)=mu_A;
err_min=1;
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
    
    out(k+1)=mu_newA;
    
    %% update best
    if out(k+1) < err_min
       err_min = out(k+1);
       best_G = G;
    end
    
    %stop condition
    if (abs(mu_newA-mu_A)/abs(mu_A)<1e-6)
        break;
    end
    mu_A=mu_newA;
    
    
end
best_mu=max(max(abs(best_G) - eye(n)));
time = toc;

%fprintf('mu_Elad=%12.8f\n',out(end));
end




function [time,best_mu] = SIDCO(A, iter)

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
err_min=1;

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
    %% update best
    if out(k+1) < err_min
       err_min = out(k+1);
       best_A = A;
    end
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
best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));

time = toc;

%fprintf('time_SIDCO=%12.8f\n mu_SIDCO=%12.8f\n iter_SIDCO=%6i\n',[time,mu,k]);

end


function [time,mu_newA] = ISPM(A, desired_coh, alpha, iter)

% Incoherence via Shifted Power Method
% Corrected version of the algorithm appeared in IEEE SPL, vol.24, Sept.2017
% "Designing Incoherent Frames With Only Matrix-ector Multiplications"
%
% Input:
%   A         - Initial frame
%   desired_coh - desired coherence, used as target (if too small, the
%                 algorithm does not converge)
%   alpha       - trade-off parameter; should be 0.5 or larger
%   iter         - number of iterations (the algorithm does not stop bases
%                 on an error criterion
% Output:
%   F           - the frame
%   out       - vector of mutual coherences at each iteration

% Note: the correction regards eq.(6): a square root should be applied on
% the right hand side

% BD 3.05.2018
tic;
% start with UNTF
F = gen_rand_untf(A);

[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
out=zeros(iter+1,1);
out(1)=mu_A;


%errv = zeros(1,Nit);
err_min = 1;

% main loop
for k = 1 : iter
  perm = randperm(n-1)+1;
  for j = perm  % a round of the atoms in random order
    d = F(:,j);
    v = F'*d;
    v(j) = 0;
    w = max(abs(v)/desired_coh, 1);
    w_max = max(w);
    eig_shift = (n/m)*alpha*w_max + (1-alpha)*sum(w)/m; % compute shift
    d = F*(w.*v) - eig_shift*d; % power method iteration
    d = d / norm(d);
    F(:,j) = d;
  end
  out(k+1) = max(max(abs(F'*F - eye(n))));
  if out(k+1) < err_min
    err_min = out(k+1);
    F_best = F;
  end
  
  %stop condition
  if (abs(out(k+1)-out(k))/abs(out(k))<1e-6)
        break;
  end
  
%   errv(k) = max(max(abs(F'*F - eye(n))));
%   if errv(k) < err_min
%     err_min = errv(k);
%     F_best = F;
%   end
  
  
  
end

% make sure that it is UNTF (most likely useless)

F = F_best;
[U,~,V] = svd(F);
F = U*V(:,1:m)'*sqrt(n/m);


mu_newA=max(max(abs(F'*F - eye(n))));
%newG=F'*F;
time = toc;
%fprintf('mu_ISPM=%12.8f\n',out(end));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = gen_rand_untf(A)

% Generate random unit norm tight frame of size mxn

% BD 10.4.2017

% generate random tight frame
%D = randn(m,n);
[m,n]=size(A);
D=A;
[Q,~] = qr(D',0);
D = Q' * sqrt(n/m);  % rows have correct norm

% force atoms to unit norm
atom_norms = sum(D.*D);
for i = 1 : n-1
  if atom_norms(i) ~= 1  % do nothing if it happens that the norm is 1
    s1 = sign(atom_norms(i)-1);
    j = i+1;
    while sign(atom_norms(j)-1) == s1  % find atom with norm on the other side of 1
      j = j+1;
    end
    % compute tangent of rotation angle
    an1 = atom_norms(i);
    an2 = atom_norms(j);
    cp = D(:,i)'*D(:,j);
    t = (cp + sign(cp)*sqrt(cp*cp - (an1-1)*(an2-1))) / (an2-1);
    % compute rotation
    c = 1 / sqrt(1+t*t);
    s = c*t;
    % new atoms and updated norm
    D(:,[i j]) = D(:,[i j]) * [c s; -s c];
    atom_norms(j) = an1+an2-1;
  end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







function [time,best_mu]=DMCM_PG(A,rho1,iter)

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
out=zeros(iter+1,1);
out(1)=mu_A;

err_min=1;

a=0.99*rho1;
for k=1:1:iter
    
    W=(A'*A-eye(n))*(1/rho1);
    v=touying_L1ball(W(:),1);
    V=reshape(v,n,n);
    Z=A-a*A*(V+V');
    A = bsxfun(@rdivide, Z, sqrt(sum(Z.^2)));
    
    
    mu_newA=max(max(abs(A'*A) - eye(n)));
    
    out(k+1)=mu_newA;
    if out(k+1) < err_min
       err_min = out(k+1);
       best_A = A;
    end
    %stop condition
    if (abs(mu_newA-mu_A)/abs(mu_A)<1e-6)
        break;
    end
    mu_A=mu_newA;
    
end
best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));

time = toc;
%fprintf('mu_DMCM-PG=%12.8f\n',out(end));
end




