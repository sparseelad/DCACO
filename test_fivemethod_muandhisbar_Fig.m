
function []=test_fivemethod_muandhisbar_Fig()

%======================================================
%This Script produces Fig. 5 and Fig. 6 in our paper.

%======================================================


%DCACO
%Elad
%SIDCO
%SIPM
%DMCM-PG

clear
clc


m = 20; % dimension of the space
n = 40; % number of vectors


app_wb=1/sqrt(m); % approximate WB
wb=sqrt((n-m)/(m*(n-1))); % WB



%% initializatio
A = randn(m, n);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));

%%improve initialization by SVD
% [U, ~, V] = svd(A);
% A = U*[eye(m) zeros(m, n-m)]*V';
% A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
% G=A'*A;
% mu_A = max(max(abs(G) - eye(n)))


G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
oldG=G;


fprintf('app_wb=%12.8f\nwb=%12.8f\nmu_A=%12.8f\n',[app_wb,wb,mu_A]);



iter =1000; % number of iterations
rho=0.2; %parameter rho

[out_DCACO,newG_DCACO]=DCACO(A,rho,iter);


dd1=0.8;
%dd21=0.2;
dd22=0.5;
dd23=0.95;

%[out_Elad1,newG_Elad1]=Elad(A,dd1,dd21,iter);
[out_Elad1,newG_Elad1]=Elad(A,dd1,dd22,iter);
[out_Elad2,newG_Elad2]=Elad(A,dd1,dd23,iter);


[out_SIDCO,newG_SIDCO] = SIDCO(A, iter);



% a_mn=3*n-m^2-2*m;
% b_mn=(m+2)*(n-m);
% lb=sqrt(a_mn/b_mn); % LB
% desired_coh=1.1*max(wb,lb);
desired_coh=1.1*wb;
alpha=0.5;

[out_ISPM,newG_ISPM] = ISPM(A, desired_coh, alpha, iter);


rho=0.5;
[out_DMCM_PG,newG_DMCM_PG] = DMCM_PG(A,rho,iter);




figure(1)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(0:1:100,mu_A*ones(101,1),'k--','LineWidth',3);hold on
plot(0:1:100,out_DCACO(1:101),'b-','LineWidth',3);hold on

plot(0:1:100,out_Elad1(1:101),'g-','LineWidth',3);hold on
plot(0:1:100,out_Elad2(1:101),'k-','LineWidth',3);hold on
plot(0:1:100,out_SIDCO(1:101),'m-','LineWidth',3);hold on
plot(0:1:100,out_ISPM(1:101),'r-','LineWidth',3);hold on
plot(0:1:100,out_DMCM_PG(1:101),'c-','LineWidth',3);hold on
plot(0:1:100,wb*ones(101,1),'r--','LineWidth',3);hold on

axis([0 100 0 1]);
xlabel('Number of iterations','FontSize',14);
ylabel('Value of \mu','FontSize',14);
legend('Initialization','DCACO','Elad (\gamma=0.5)','Elad (\gamma=0.95)','SIDCO','ISPM','DMCM-PG','Welch bound')


figure(2)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);
plot(0:1:iter,mu_A*ones(iter+1,1),'k--','LineWidth',3);hold on
plot(0:1:iter,out_DCACO,'b-','LineWidth',3);hold on

plot(0:1:iter,out_Elad1,'g-','LineWidth',3);hold on
plot(0:1:iter,out_Elad2,'k-','LineWidth',3);hold on
plot(0:1:iter,out_SIDCO,'m-','LineWidth',3);hold on
plot(0:1:iter,out_ISPM,'r-','LineWidth',3);hold on
plot(0:1:iter,out_DMCM_PG,'c-','LineWidth',3);hold on
plot(0:1:iter,wb*ones(iter+1,1),'r--','LineWidth',3);hold on

axis([0 iter 0 1]);
xlabel('Number of iterations','FontSize',14);
ylabel('Value of \mu','FontSize',14);
legend('Initialization','DCACO','Elad (\gamma=0.5)','Elad (\gamma=0.95)','SIDCO','ISPM','DMCM-PG','Welch bound')


g1=sort(abs(oldG(:)));
oldg1=g1(1:1:end-n);
[a1,b1]=hist(oldg1,0:0.01:1);

gDCACO=sort(abs(newG_DCACO(:)));
newgDCACO=gDCACO(1:1:end-n);
[a2,b2]=hist(newgDCACO,0:0.01:1);



gElad1=sort(abs(newG_Elad1(:)));
newgElad1=gElad1(1:1:end-n);
[a3,b3]=hist(newgElad1,0:0.01:1);

gElad2=sort(abs(newG_Elad2(:)));
newgElad2=gElad2(1:1:end-n);
[a4,b4]=hist(newgElad2,0:0.01:1);


gSIDCO=sort(abs(newG_SIDCO(:)));
newgSIDCO=gSIDCO(1:1:end-n);
[a5,b5]=hist(newgSIDCO,0:0.01:1);


gISPM=sort(abs(newG_ISPM(:)));
newgISPM=gISPM(1:1:end-n);
[a6,b6]=hist(newgISPM,0:0.01:1);

gDMCM_PG=sort(abs(newG_DMCM_PG(:)));
newgDMCM_PG=gDMCM_PG(1:1:end-n);
[a7,b7]=hist(newgDMCM_PG,0:0.01:1);


figure(3)
subplot(7,1,1)
bar(b1,a1,'b')
legend('Initialization')
axis([0 1 -inf inf]);
subplot(7,1,2)
bar(b2,a2,'b')
legend('DCACO')
axis([0 1 -inf inf]);



subplot(7,1,3)
bar(b3,a3,'b')
legend('Elad (\gamma=0.5)')
axis([0 1 -inf inf]);

subplot(7,1,4)
bar(b4,a4,'b')
legend('Elad (\gamma=0.95)')
axis([0 1 -inf inf]);

subplot(7,1,5)
bar(b5,a5,'b')
legend('SIDCO')
axis([0 1 -inf inf]);


subplot(7,1,6)
bar(b6,a6,'b')
legend('ISPM')
axis([0 1 -inf inf]);


subplot(7,1,7)
bar(b7,a7,'b')
legend('DMCM-PG')
axis([0 1 -inf inf]);




end








function [out,best_G]=DCACO(A,rho,iter)

% =====================================================
% also see the Script of DCACO in our software package.
% ======================================================


[m,n]=size(A);
A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
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
    
    
end

best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));
fprintf('mu_DCACO=%12.8f\n',best_mu);
end


function [out,best_G]=Elad(A,dd1,dd2,iter)

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
    
    
    %% shrink the high inner products
    % pos=find(abs(G(:))>=Threshold & abs(G(:)-1)>1e-6);
    % G(pos)=G(pos)*dd2;
    
    gg=sort(abs(G(:)));
    pos=find(abs(G(:))>gg(round(dd1*(n^2-n))) & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    
    
    %% reduce the rank back to m
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
    
end

best_mu=max(max(abs(best_G) - eye(n)));
fprintf('mu_Elad=%12.8f\n',best_mu);
end








function [out,best_G] = SIDCO(A, iter)

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
        
        
        newA=PC_sub(Awork',target,sqrt(1-T^2));
        
        newA = newA/norm(newA);
        A(:, i) = newA;
        
    end
    %% new frame, new mutual coherence
    mu = max(max(abs(A'*A) - eye(n)));
    out(k+1)=mu;
    %% update best
    if out(k+1) < err_min
       err_min = out(k+1);
       best_A = A;
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
%A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));
fprintf('mu_SIDCO=%12.8f\n',best_mu);
end




function [out,best_G]=DMCM_PG(A,rho,iter)

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





[~,n]=size(A);

A = bsxfun(@rdivide, A, sqrt(sum(A.^2)));
G=A'*A;
mu_A=max(max(abs(G) - eye(n)));
out=zeros(iter+1,1);
out(1)=mu_A;

err_min=1;


a=0.99*rho;
for k=1:1:iter
    
    W=(A'*A-eye(n))*(1/rho);
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
    
end
best_G=best_A'*best_A;
best_mu=max(max(abs(best_G) - eye(n)));
fprintf('mu_DMCM-PG=%12.8f\n',best_mu);
end







function [out,best_G] = ISPM(A, desired_coh, alpha, iter)

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



best_G=F'*F;
best_mu=max(max(abs(best_G) - eye(n)));
fprintf('mu_ISPM=%12.8f\n',best_mu);
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
