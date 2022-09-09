function P=DesignProjection_Elad_D(D,D_typ,n,Iter,dd1,dd2,Init)

%======================================================================
% This Script is directly from the software package: "SparseLab
% \SparseLab2.1-Core\CompSense" which can be downloaded free from the 
% website: https://elad.cs.technion.ac.il/software/
% 
% More details on the main idea of this Script, please see the paper: 
% "Optimized Projections for Compressed Sensing," 
% IEEE Trans. Signal Process., vol. 55, no. 12, pp. 5695-5702, 2007.

%======================================================================



%========================================
% In this program we iteratively build a projection matrix that when
% applied to the dictionary D, it leads to better coherence.
%
% Inputs:    D - The dictionary
%               n - Number of projections (rows in P)
%               Iter - Number of iterations
%               dd1 - Relative number of the G entries to shrink
%               dd2 - Shrink factor to use (1/dd2 is used to expand)
%               Init - An initial projection matrix
%
% Outputs:  P - The projection matrix
%
% Example: N=100; L=200;
%               D=randn(N,L);
%               D=D*diag(1./sqrt(diag(D'*D)));
%               n=20; Iter=100; dd1=0.7; dd2=0.95;
%               P=randn(n,N);
%               P=DesignProjection_Elad_D(D,n,Iter,dd1,dd2,P);
%========================================

[N,L]=size(D);
%disp(['The best achievable \mu is ',num2str(sqrt((L-n)/(n*(L-1))))]);
% pause;

if nargin==5,
    P=randn(n,N); % initialization
    % P=P*diag(1./sqrt(diag(P'*P))); % normalize columns
else
    P=Init;
end
G=D'*P'*P*D; % compute the Gram matrix of the projected dictionary
G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G))); % nromalize columns
% gg=sort(abs(G(:)));
% pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
% RefVal=mean(abs(G(pos)));

for k=1:1:Iter,
    
    % shrink the high inner products
    gg=sort(abs(G(:)));
    pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    G(pos)=G(pos)*dd2;
    
    % expand the near zero products
    % pos=find(abs(G(:))<gg(round((1-dd1)*(L^2-L))));
    % G(pos)=G(pos)/dd2;
    
    % reduce the rank back to n
    method=2;
    if method==1,
        [U,S,V]=svds(G,n);
        S=S(1:n,1:n);
        U=U(:,1:n);
    elseif method==2,
        U=randn(L,n);
        for jjj=1:1:10, U=G*U; U=orth(U); end;
        S=diag(diag(U'*G*U));
    end;
    G=U*S*U';
    PD=S.^0.5*U';
    
    if D_typ==1  
       P=PD;
    else
       P=PD*pinv(D);
    end
    
    %P=PD*pinv(D);
    
    % Normalize the columns
    G=D'*P'*P*D;
    G=diag(1./sqrt(diag(G)))*G*diag(1./sqrt(diag(G)));
    
    % Show status
    %     gg=sort(abs(G(:)));
    %     pos=find(abs(G(:))>gg(round(dd1*(L^2-L))) & abs(G(:)-1)>1e-6);
    %     if RefVal>mean(abs(G(pos))),
    %         RefVal=mean(abs(G(pos)));
    %     else,
    %         return;
    %     end;
    %     fprintf(1,'%6i  %12.8f  %12.8f \n',[k,mean(abs(G(pos))),max(abs(G(pos)))]);
    
end

end
