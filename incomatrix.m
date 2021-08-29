function A=incomatrix(X,m)
%==================================================
%This Script computes an m x n matrix from a given 
%symmetric matrix X by using svd(X).

%===================================================

% [U,S]=eig(X);% use eig
[U,S,V]=svd(X);%Instead of using eig, we use svd to guarantee that the diagonal elements are nonnegative.
S=S(1:m,:);
DT=sqrt(S);
A=DT*U';

end