

%====================================================
%This Script produces Fig. 8 in our paper.

%====================================================



clear;
clc;


% Setting Parameters
K=120; % number of atoms in the dictionary
n=80; % dimension of the signals
q=100000; % number of examples
%qT=300; % number of error to accumulate
m=25; % number of projections
% Creation of the dictionary
D_typ=0;
D=randn(n,K);
D=D*diag(1./sqrt(diag(D'*D)));

% Generating the compressed measurements
P=randn(m,n);
for j=1:1:m,
    P(j,:)=P(j,:)/norm(P(j,:));
end



iter=1000; dd1=0.8; dd2=0.95;
rho=0.2;
bate=2;

wb=sqrt((n-m)/(m*(n-1))); % WB

desired_coh=1.1*wb;
alpha=0.5;



% optimizing the projection  by DCACO, Elad and DMCM-AM


DCACO_Pnew=DesigenProjection_DCACO_D(P,D,D_typ,rho,iter);

Elad_Pnew=DesignProjection_Elad_D(D,D_typ,m,iter,dd1,dd2,P);

SIDCO_Pnew=DesigenProjection_SIDCO_D(P,D,D_typ,iter);

ISPM_Pnew=DesigenProjection_ISPM_D(P,D,D_typ,desired_coh, alpha, iter);

DMCM_AM_Pnew=DesigenProjection_DMCM_AM_D(P,D,D_typ,rho,bate,iter);

ss=7;
RESULT=zeros(ss,6);
for s=1:1:ss, % cardinality of the sparse representation
    
    % Creation of test signals having sparse representations
    A=zeros(K,q); A=sparse(A);
    h=waitbar(0,'Building the sparse represnetations ...');
    for j=1:1:q,
        waitbar(j/q);
        pos=randperm(K);
        pos=pos(1:s);
        A(pos,j)=randn(s,1);
    end
    close(h);
    X=D*A;
    
    
    %  decoding using OMP
    Ahat=A*0;
    Dhat=P*D; Y=P*X;
    out=zeros(q,1);
    for j=1:1:q,
        Ahat(:,j)=OMP(Dhat,Y(:,j),s);
        out(j)=mean((Ahat(:,j)-A(:,j)).^2)/mean(A(:,j).^2)>1e-3;
        %         if sum(out(1:j))>=qT,
        %             break;
        %         end;
    end
    RESULT(s,1)=mean(out(1:j));
    
    
    
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
    RESULT(s,2)=mean(out(1:j));
    
    
    
    
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
    RESULT(s,3)=mean(out(1:j));
    
  
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
    RESULT(s,4)=mean(out(1:j));  
    
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
    RESULT(s,5)=mean(out(1:j));  
    
    
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
    RESULT(s,6)=mean(out(1:j));
    
    fprintf('s=%d\n',s);
end



figure(1)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);

Y1=RESULT(:,1); Y2=RESULT(:,2);Y3=RESULT(:,3);Y4=RESULT(:,4);Y5=RESULT(:,5);Y6=RESULT(:,6);
plot(1:1:ss,Y1,'-*k','LineWidth',2,'MarkerSize',10); hold on; plot(1:1:ss,Y2,'-*b','LineWidth',2,'MarkerSize',10);
hold on; plot(1:1:ss,Y3,'-*g','LineWidth',2,'MarkerSize',10); hold on; plot(1:1:ss,Y4,'-*m','LineWidth',2,'MarkerSize',10);
hold on; plot(1:1:ss,Y5,'-*r','LineWidth',2,'MarkerSize',10);hold on; plot(1:1:ss,Y6,'-*c','LineWidth',2,'MarkerSize',10);

legend('Random','DCACO','Elad','SIDCO','ISPM','DMCM-AM','Location','NorthWest')
xlabel('Sparsity','FontSize',14);
ylabel('Relative # of errors','FontSize',14);
%grid on;