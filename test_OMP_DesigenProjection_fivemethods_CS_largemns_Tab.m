

%==========================================================================
%This Script produces results in TABLE 5, TABLE 6 and TABLE 7 in our paper.

%==========================================================================




clear;
clc;


% Setting Parameters
%K=120; % number of atoms in the dictionary
n=2000; % dimension of the signals
q=50; % number of examples
%qT=300; % number of error to accumulate
%s=4; % cardinality of the sparse representation
%m=16:2:40; % number of projections

D_typ=1; %type of dictionary

RESULT=zeros(length(n),6); %Initialization result
for jjj=1:1:length(n), % using different number of projections
    
    
    K=n(jjj);% number of atoms in the dictionary
    m=n(jjj)/10;% number of projections
    s=m/5;% cardinality of the sparse representation
    
    % Creation of the dictionary
    D=eye(n(jjj));
    
    
    % Creation of test signals having sparse representations
    A=zeros(K,q); A=sparse(A);
    h=waitbar(0,'Building the sparse representations ...');
    for j=1:1:q,
        waitbar(j/q);
        pos=randperm(K);
        pos=pos(1:s);
        A(pos,j)=randn(s,1);
    end
    close(h);
    X=D*A;
    
    
    % Generating the compressed measurements
    
    P=randn(m,n(jjj));
    for j=1:1:m,
        P(j,:)=P(j,:)/norm(P(j,:));
    end
    Y=P*X; % the measurements
    
    
    %  decoding using OMP
    Ahat=A*0;
    Dhat=P*D;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,1)=mean(out(1:j));
    
    
    
    % optimizing the projection  by DCACO
    iter=500; rho=0.2;
    DCACO_Pnew=DesigenProjection_DCACO_D(P,D,D_typ,rho,iter);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=DCACO_Pnew*D;
    Y=DCACO_Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,2)=mean(out(1:j));
    
    
    
    % optimizing the projection by Elad
    iter=500; dd1=0.8; dd2=0.95;
    Elad_Pnew=DesignProjection_Elad_D(D,D_typ,m,iter,dd1,dd2,P);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=Elad_Pnew*D;
    Y=Elad_Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,3)=mean(out(1:j));
    
    % optimizing the projection  by SIDCO
    iter=100; 
    
    SIDCO_Pnew=DesigenProjection_SIDCO_D(P,D,D_typ,iter);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=SIDCO_Pnew*D;
    Y=SIDCO_Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,4)=mean(out(1:j));
    
 
    
    
    % optimizing the projection  by ISPM
    iter=1000; 
    wb=sqrt((n(jjj)-m)/(m*(n(jjj)-1))); % WB
    desired_coh=1.1*wb;
    alpha=0.5;
    
    ISPM_Pnew=DesigenProjection_ISPM_D(P,D,D_typ,desired_coh, alpha, iter);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=ISPM_Pnew*D;
    Y=ISPM_Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,5)=mean(out(1:j));
    
    
    
    % optimizing the projection  by DMCM-AM
    iter=500; rho=0.2;bate=2;
    DMCM_AM_Pnew=DesigenProjection_DMCM_AM_D(P,D,D_typ,rho,bate,iter);
    
    % the greedy new results
    Ahat=A*0;
    Dhat=DMCM_AM_Pnew*D;
    Y=DMCM_AM_Pnew*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(jjj,6)=mean(out(1:j));
    
    fprintf('jjj=%d\n',jjj);
    
end
    
    
disp(RESULT)
    
save('RESULT.mat','RESULT')  