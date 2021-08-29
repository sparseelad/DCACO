function mu=compu_mu(A)
%=======================================================
%This Script computes the mutual coherence of a matrix A.

%=======================================================

n=size(A,2);



% for i=1:1:n
%     A(:,i)=A(:,i)/norm(A(:,i));
% end
%G=A'*A-eye(n);
%mu=norm(G(:),Inf);

A=bsxfun(@rdivide, A, sqrt(sum(A.^2))); % Normalize the columns
mu=max(max(abs(A'*A) - eye(n)));
end