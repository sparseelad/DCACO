function y=sub_diff(X,m)

%=======================================================================
%This Script chooses an element y in the subdifferential set of the 
%matrix convex function \Phi at the point x.

%=======================================================================


[Q,D]=eig(X);
d=diag(D);
dd=sort(d,'descend');
ddm=dd(m);
s=find(d>ddm);

if size(s,1)==0
    t=find(d==ddm);
    tm=t(1:m);
    d_new=max(d(tm),0);
else
    
    if ddm<=0
        tm=s;
        d_new=max(d(tm),0);
    else
        t=find(d==ddm);
        tt=[s;t];
        tm=tt(1:m);
        d_new=d(tm);
    end
    
end
Q=Q(:,tm);
y=Q*diag(d_new)*Q';

end

