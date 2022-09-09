

%====================================================
%This Script produces Fig. 7 in our paper.

%====================================================



clear;
clc;


% Setting Parameters
K=120; % number of atoms in the dictionary
n=80; % dimension of the signals
q=100000; % number of examples
%qT=300; % number of error to accumulate
s=4; % cardinality of the sparse representation
m=16:2:40; % number of projections

% Creation of the dictionary
D_typ=0;
D=randn(n,K);
D=D*diag(1./sqrt(diag(D'*D)));

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

RESULT=zeros(length(m),6);
for jjj=1:1:length(m), % using different number of projections
    
    % Generating the compressed measurements
    P=randn(m(jjj),n);
    for j=1:1:m(jjj),
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
    iter=1000; rho=0.2;
    
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
    iter=1000; dd1=0.8; dd2=0.95;
    Elad_Pnew=DesignProjection_Elad_D(D,D_typ,m(jjj),iter,dd1,dd2,P);
    
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
    iter=1000; 
    
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
    
 
    
    
    % optimizing the projection  by SIDCO
    iter=1000; 
    wb=sqrt((n-m(jjj))/(m(jjj)*(n-1))); % WB
%     a_mn=3*n-(m(jjj))^2-2*m(jjj);
%     b_mn=(m(jjj)+2)*(n-m(jjj));
%     lb=sqrt(a_mn/b_mn); % LB
%     desired_coh=1.1*max(wb,lb);
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
    iter=1000; rho=0.2;bate=2;
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


figure(1)
set(gca,'FontSize',14,'position',[0.35812499999999997,0.25264270613107837,0.3080208333333333,0.625792811839324]);

Y1=RESULT(:,1); Y2=RESULT(:,2);Y3=RESULT(:,3);Y4=RESULT(:,4);Y5=RESULT(:,5);Y6=RESULT(:,6);
plot(m,Y1,'-*k','LineWidth',2,'MarkerSize',10); hold on; plot(m,Y2,'-*b','LineWidth',2,'MarkerSize',10);
hold on; plot(m,Y3,'-*g','LineWidth',2,'MarkerSize',10); hold on; plot(m,Y4,'-*m','LineWidth',2,'MarkerSize',10);
hold on; plot(m,Y5,'-*r','LineWidth',2,'MarkerSize',10);hold on; plot(m,Y6,'-*c','LineWidth',2,'MarkerSize',10);

legend('Random','DCACO','Elad','SIDCO','ISPM','DMCM-AM')
xlabel('Number of measurements','FontSize',14);
ylabel('Relative # of errors','FontSize',14);
%grid on;