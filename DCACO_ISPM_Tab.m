
%==================================================================
%This Script produces results in TABLE 4 in our paper.

%==================================================================



function []=DCACO_ISPM_Tab()

clear
clc


m = 20; % dimension of the space
n = 40; % number of vectors
% m=10;
% n=100;
% m=20;
% n=500;
%app_wb=1/sqrt(m); % approximate WB
wb=sqrt((n-m)/(m*(n-1))); % WB

% a_mn=3*n-m^2-2*m;
% b_mn=(m+2)*(n-m);
% lb=sqrt(a_mn/b_mn); % LB
% desired_coh=1.1*max(wb,lb);
desired_coh=1.1*wb;


alpha=0.5;

rho=0.2;

iter=1000;
iter_ISPM=1000;
%% initializatio
A = randn(m, n);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));

[time1,mu1,F,k1]=ISPM(A, desired_coh, alpha, iter_ISPM);

[time2,mu2,k2]=DCACO(A,rho,iter);


[time3,mu3,k3]=DCACO(F,rho,iter);

fprintf('DCACO_time=%12.8f\nISPM_time=%12.8f\nDCACObyISPM_time=%12.8f\n',[time1,time2,time3]);
fprintf('DCACO_iternum=%6i\nISPM_iternum=%6i\nDCACObyISPM_iternum=%6i\n',[k1,k2,k3]);
fprintf('wb=%12.8f\nmu_A=%12.8f\nDCACO_mu=%12.8f\nISPM_mu=%12.8f\nDCACObyISPM_mu=%12.8f\n',[wb,mu_A,mu1,mu2,mu3]);
end


function [time,best_mu,k]=DCACO(A,rho,iter)

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
err_min = 1;
for k=1:1:iter
    Y=sub_diff(G,m);
    W=Y-eye(n);
    w=smatrix_to_vec(W);
    u=touying_L1ball(w,1/(2*rho));
    z=w-u;
    G=vec_to_smatrix(z,n);
    
    newA=incomatrix(G,m);
    %mu_newA=compu_mu(newA);
    newA = bsxfun(@rdivide, newA, sqrt(sum(newA.^2)));
    newG=newA'*newA;
    mu_newA=max(max(abs(newG) - eye(n)));
    
    out(k+1)=mu_newA;
    
    %% update best
    if out(k+1) < err_min
       err_min = out(k+1);
       best_A = newA;
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

%fprintf('time_DCACO=%12.8f\n mu_DCACO=%12.8f\n iter_DCACO=%6i\n',[time,mu_newA,k]);

end







function [time,mu_newA,F,k] = ISPM(A, desired_coh, alpha, iter)

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

[m,n]=size(F);

errv=zeros(iter,1);


err_min = 1;
w_maxmax=1;
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
    
    w_maxmax = max(w_maxmax,sqrt(w_max));% the maxmum weight in inner loop
  end
  errv(k) = max(max(abs(F'*F - eye(n))));
  if errv(k) < err_min
    err_min = errv(k);
    F_best = F;
  end
  
  %stop condition
  if (w_maxmax<=1.0001)
        break;
  end
  
   
end

% make sure that it is UNTF (most likely useless)

F = F_best;
[U,~,V] = svd(F);
F = U*V(:,1:m)'*sqrt(n/m);


mu_newA=max(max(abs(F'*F - eye(n))));

time = toc;

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