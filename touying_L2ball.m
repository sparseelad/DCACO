function u=touying_L2ball(w,y,r)

% ==================================================================
% 
% 
% compute the projection operator P_C onto the closed convex set
%                     
%                       C={x:||x-y||<=r}
%                         =y+{x:||x||_2<=r}.
% 
% ===================================================================                      
                      
                      
a=r/max(norm(w-y),r);
u=y+a*(w-y);


end





