function newP=DesigenProjection_ISPM_D(P,D,D_typ,desired_coh, alpha, iter)

% ===========================================================================
% This Script produces an optimized projection matrix by applying ISPM when
% given a dictionary D.

%
% First, applying ISPM with an
% initial effective dictionary A = PD to produce an incoherent
% effective dictionary A_ISPM.
% Second, solving a least squares
% problem, we obtain an optimized projection matrix P_ISPM
% satisfying A_ISPM is approximately equal to P_ISPM.
% 
% ===========================================================================

%%ISPM
A=P*D;

% start with UNTF
F = gen_rand_untf(A);

[m,n]=size(F);
errv = zeros(1,iter);
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
  errv(k) = max(max(abs(F'*F - eye(n))));
  if errv(k) < err_min
    err_min = errv(k);
    F_best = F;
  end
end
% make sure that it is UNTF (most likely useless)
F = F_best;
[U,~,V] = svd(F);
F = U*V(:,1:m)'*sqrt(n/m);

%%solving a least squares problem
if D_typ==1  
newP=F;
else
newP=F*pinv(D);
end

%newP=F*pinv(D);


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