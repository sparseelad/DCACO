function newP=DesigenProjection_SIDCO_D(P,D,D_typ,iter)

% ===========================================================================
% This Script produces an optimized projection matrix by applying ISPM when
% given a dictionary D.

%
% First, applying SIDCO with an
% initial effective dictionary A = PD to produce an incoherent
% effective dictionary A_SIDCO.
% Second, solving a least squares
% problem, we obtain an optimized projection matrix P_SIDCO
% satisfying A_SIDCO is approximately equal to P_SIDCO.
% 
% ===========================================================================

%%SIDCO
A=P*D;



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

%%solving a least squares problem
%newP=A*pinv(D);
if D_typ==1   
   newP=best_A;
else
   newP=best_A*pinv(D);
end


end