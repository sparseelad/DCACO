function w=smatrix_to_vec_simple(W)

%==================================================
%This Script constructs the vector w from a given 
%symmetric matrix W. 

%It corresponds to the function 'smatrix_to_vec'.

%===================================================



n=size(W,1);
w=[];

for i=1:1:(n-1)
    
    w1=[w,W(i,(i+1):n)];
    w=w1;
    
end
w=w';
end