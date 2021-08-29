function w=smatrix_to_vec(W)

%==================================================
%This Script constructs the vector w from a given 
%symmetric matrix W.

%===================================================


n=size(W,1);
q=n*(n-1);
q=q/2;
w=zeros(q,1);
for i=1:(n-1)
    
    a=(i-1)*(2*n-i);
    a=1+a/2;
    b=a+n-(i+1);
    w(a:b)=W(i,(i+1):n)';
    
end
end