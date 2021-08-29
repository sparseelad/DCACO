function u=touying_L1ball(w,r)

% ==================================================================
% 
% 
% compute the projection operator P_D onto the closed convex set
%                     
%                       C={z:||z||_1<=r}.
%                         
% 
% ===================================================================


n=length(w);
t=norm(w,1);
if t<=r
u=w;
else
v=abs(w);
wd=sort(v,'descend');
i=1;
tau=wd(1)-r;

while (wd(i)-tau+eps)>0
    J=i; %Record the value of the previous i.
    Jtau=tau; %Record the value of the previous tau.
    if (i==n)
        break;
    end
    i=i+1; %Update the value of i.
    tau=(sum(wd(1:i))-r)/i; %Update the value of tau.
end

u=sign(w).*max(v-Jtau,0);

end

end