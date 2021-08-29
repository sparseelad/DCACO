function a=OMP(D,x,L)


%======================================================================
% This Script implements Orthonormal Matching Pursuit with L non-zeros.

% This Script is directly from Fig. 3.1 in Chapter 3 in 
% M. Elad's book: "Sparse and Redundant Representations: From Theory to
% Applications in Signal and Image Processing," Springer, New York, 2010.
% 

%======================================================================


[n,K]=size(D);
a=[];
residual=x;
indx=zeros(L,1);
for j=1:1:L,
    proj=D'*residual;
    pos=find(abs(proj)==max(abs(proj)));
    pos=pos(1);
    indx(j)=pos;
    a=pinv(D(:,indx(1:j)))*x;
    residual=x-D(:,indx(1:j))*a;
end;
temp=zeros(K,1); 
temp(indx)=a;
a=sparse(temp);

end