%%%%%%%%%%%%%%
% THRESH_P.M %
%%%%%%%%%%%%%%
% This is a Matlab program file.
% It replicates the estimation, testing and graphs reported in
% "Threshold Effects in Non-Dynamic Panels:
% Estimation, Testing and Inference"

% For questions, please contact

% Bruce E. Hansen
% Department of Economics
% Social Science Building
% University of Wisconsin
% Madison, WI 53706-1393
% bhansen@ssc.wisc.edu
% http://www.ssc.wisc.edu/~bhansen/


% This program file loads the GAUSS dataset "invest.fmt".
% It creates the output file "thresh.out"
function thresh_p;
global nt;
global t;
global n;
global max_lag;
global thresh;
global cf;
global qn;
global qq1;
global vgraph_;
global tt;
global yt;
global xt;
global ct;
global ty;
global cc;
global k;
global tt;
global qn1;

load invest.txt;
t = 15;
nt = length(invest(:,1));
n = nt/t;

i = invest(:,1);     % investment/assets
q = invest(:,2);     % Tobin's Q    
c = invest(:,3);     % cash-flow/assets                              
d = invest(:,4);     % debt/assets

qn = 20;%400;            % number of quantiles to examine                %
conf_lev = .95;      % confidence level for threshold                %
vgraph_ = 1;         % set to 1 to graph likelihood ratios           %
boot_1_ = 50;%300;       % # of replications, 0 for no bootstrap, single (300) %
boot_2_ = 20;%300;       % # of replications, 0 for no bootstrap, double (300) %
boot_3_ = 10;%300;       % # of replications, 0 for no bootstrap, triple (300) %
trim_1_ = .01;       % percentage to trim before search, single      %
trim_2_ = .01;       % percentage to trim before search, double      %
trim_3_ = .05;       % percentage to trim before search, triple      %

max_lag = 1;
tt = t-max_lag;
ty = n*(t-max_lag-1); % AT: Why (t-max_lag-1) ?

y  = lag_v(i,0);  yt = tr(y);
cf = lag_v(c,1);  ct = tr(cf);
q1 = lag_v(q,1);
d1 = lag_v(d,1);      % set to threshold variable %

x=[q1,q1.^2,q1.^3,d1,q1.*d1];
k=length(x(1,:));
xt=zeros(length(yt(:,1)),k);
j=1;
while j<=k
    xt(:,j)=tr(x(:,j));
    j=j+1;
end;

thresh=d1;
dd=unique(thresh);
qnt1=qn*trim_1_;
sq=trim_1_:(1/qn):(trim_1_+(1/qn)*(qn-2*qnt1));
qq1=dd(floor(sq*length(dd(:,1))));
qn1=length(qq1(:,1));
cc=-2*log(1-sqrt(conf_lev));
out=fopen('out.txt','wt');
fprintf(out,'Number of Firms:                 %u\n',n);
fprintf(out,'Number of Years used:            %u\n',tt);
fprintf(out,'Total Observations:              %u\n',ty);
fprintf(out,'Number of Quantiles:             %u\n',qn);
fprintf(out,'Confidence Level:                %f\n',conf_lev);
fprintf(out,'\n');
fprintf(out,'\n');
fprintf(out,'*********************************************\n');
fprintf(out,'\n');
fprintf(out,'\n');

sse0 = sse_calc(yt,[xt,ct]);
fprintf(out,'Zero Threshold Model\n');
fprintf(out,'Sum of Squared Errors:           %f\n',sse0);
fprintf(out,'\n');
fprintf(out,'\n');
fprintf(out,'*********************************************\n');
fprintf(out,'\n');
fprintf(out,'\n');

fprintf(out,'Single Threshold Model\n');
fprintf(out,'\n');
fclose(out);
rhat1=model(0,trim_1_,boot_1_,0);
out=fopen('out.txt','a+');
fprintf(out,'*********************************************\n');
fprintf(out,'\n');
fprintf(out,'\n');

fprintf(out,'Double Threshold Model\n');
fprintf(out,'Trimming Percentage:             %f\n',trim_2_);
fprintf(out,'\n');
fprintf(out,'First Iteration\n');
fclose(out);
rhat2=model(rhat1,trim_2_,boot_2_,2);
out=fopen('out.txt','a+');
fprintf(out,'Second Iteration\n');
fclose(out);
rhat1=model(rhat2,trim_2_,0,1);
out=fopen('out.txt','a+');
fprintf(out,'\n');
fprintf(out,'\n');
fprintf(out,'*********************************************\n');
fprintf(out,'\n');
fprintf(out,'\n');

fprintf(out,'Triple Threshold Model\n');
fprintf(out,'Trimming Percentage:             %f\n',trim_3_);
fprintf(out,'\n');
fclose(out);
rhat3=model([rhat1;rhat2],trim_3_,boot_3_,3);
out=fopen('out.txt','a+');
fprintf(out,'\n');
fprintf(out,'\n');
fprintf(out,'*********************************************\n');
fprintf(out,'\n');
fprintf(out,'\n');
fclose(out);

% Functions %
function r=tr(y);
global n;
global tt;

yf=reshape(y',tt,n)';
yfm=yf-mean(yf')'*ones(1,length(yf(1,:)));
yfm=yfm(:,1:tt-1)'; % AT: Why (:,1:tt-1) ?
ind=0;
for i=1:length(yfm(1,:))
    for j=1:length(yfm(:,1))
        if ind==0
            r=yfm(j,i);
            ind=1;
        else
            r=[r;yfm(j,i)];
        end;
    end;
end;

function ss=sse_calc(y,x)
e=y-x*(y'/x')';
ss=e'*e;

function r=lag_v(x,lagn);
global max_lag;
global t;
global n;
y=reshape(x',t,n)';
y=y(:,(1+max_lag-lagn):(t-lagn))';
ind=0;
for i=1:length(y(1,:))
    for j=1:length(y(:,1))
        if ind==0
            r=y(j,i);
            ind=1;
        else
            r=[r;y(j,i)];
        end;
    end;
end;

function sse=thr_sse(y,q,r);
global thresh;
global cf;
global xt;
global ct;
n=length(q(:,1));
sse=zeros(n,1);
qi=1;
while qi<=n
    if r==0
        rr=q(qi);
    else
        rr=[r;q(qi)];
    end;
    rr=sortrows(rr,1);
    xx=[xt,ct];
    j=1;    
    while j<=length(rr(:,1));
        d=(thresh<rr(j));
        xx=[xx,tr(cf.*d)];
        j=j+1;
    end;
    sse(qi)=sse_calc(y,xx);
    qi=qi+1;
end;

function [rsse,rqq]=r_est(y,r,trim_);
global qn;
global qn1;
global qq1;
if max(r)'==0
    qq=qq1;
    rr=0;
else
    rr=sortrows(r,1);
    i=(1:qn1)';    
    nn=sum((qq1*ones(1,length(rr)))<(ones(length(qq1(:,1)),1)*(rr')));
    qnt=qn*trim_;
    ii=(((i*ones(1,length(nn)))<=(ones(length(i(:,1)),1)*(nn+qnt)))-((i*ones(1,length(nn)))<=(ones(length(i(:,1)),1)*(nn-qnt))))*ones(length(rr(:,1)),1);
    ind=0;
    for j=1:length(ii)
        if ii(j)==0
            if ind==0
                qq=qq1(j);
                ind=1;
            else
                qq=[qq;qq1(j)];
            end;
        end;
    end;
end;
sse=thr_sse(y,qq,rr);
[temp,rihat]=min(sse);
clear temp;
rsse=sse(rihat);
rqq=qq(rihat);

function rhat=model(r,trim_,rep,it);
global qq1;
global vgraph_;
global yt;
global ty;
global cc;
global xt;
global thresh;
global cf;
global n;
global k;
global ct;
global tt;
global qn1;
global qn;

if max(r)==0
    qq=qq1;
    rr=0;
else
    rr=sortrows(r,1);
    i=(1:qn1)';
    ind=0;
    for ij=1:length(rr)
        if ind==0
           nn=sum(qq1<rr(ij));
           ind=1;
        else
            nn=[nn,sum(qq1<rr(ij))];
        end;
    end;    
    qnt=qn*trim_;
    inn1=(i*ones(1,length(nn))<=(ones(length(i),1)*(nn+qnt)));
    inn2=(i*ones(1,length(nn))<=(ones(length(i),1)*(nn-qnt)));
    ii=(inn1-inn2)*ones(length(rr(:,1)),1);
    ind=0;
    for j=1:length(ii)
        if ii(j)==0
            if ind==0
                qq=qq1(j);
                ind=1;
            else
                qq=[qq;qq1(j)];
            end;
        end;
    end;
end;

sse=thr_sse(yt,qq,rr);
[temp,rihat]=min(sse);
clear temp;
rhat=qq(rihat);
sse1=sse(rihat);
lr=(sse/sse1-1)*ty;
ind=0;
temp=(lr<cc);
for j=1:length(qq(:,1))
    if temp(j)==1
        if ind==0
            rhats=qq(j);
            ind=1;
        else
            rhats=[rhats;qq(j)];
        end;
    end;
end;
if vgraph_==1
    figure;
    plot(qq,lr,qq,ones(length(qq),1)*cc);
    if it==0
        title('Figure 1 LConfidence Interval Construction in Single Threshold Model');
        xlabel('Threshold Parameter');
    end;
    if it==1
        title('Figure 3 LConfidence Interval Construction in Double Threshold Model');
        xlabel('First Threshold Parameter');
    end;
    if it==2
        title('Figure 2 LConfidence Interval Construction in Double Threshold Model');
        xlabel('Second Threshold Parameter');  
    end;
    if it==3
        title('Confidence Interval Construction in Triple Threshold Model');
        xlabel('Third Threshold Parameter');  
    end;
    ylabel('Likelihood Ratio');                
end;

out=fopen('out.txt','a+');     
if abs(max(r)')>0
    fprintf(out,'Fixed Thresholds:   \n');
    for i=1:length(rr) 
        fprintf(out,'%f   ',rr(i));
    end;
    fprintf(out,'\n');
    rrr=sortrows([rr;rhat],1);
else
    rrr=rhat;
end;

fprintf(out,'Threshold Estimate:              %f\n',rhat);
fprintf(out,'Confidence Region:               %f   %f\n',min(rhats),max(rhats));
fprintf(out,'Sum of Squared Errors:           %f\n',sse1);
fprintf(out,'Trimming Percentage:             %f\n',trim_);
fprintf(out,'\n');
fprintf(out,'\n');

nr=length(rrr(:,1));
xx=xt;
dd=zeros(length(thresh(:,1)),nr);
j=1;
while j<=nr
    dd(:,j)=(thresh<rrr(j));
    d=dd(:,j);
    if j>1;
        d=d-dd(:,j-1);
    end;
    xx=[xx,tr(cf.*d)];
    j=j+1;
end;
d=1-dd(:,nr);
xx=[xx,tr(cf.*d)];
xxi=inv(xx'*xx);
beta=xxi*(xx'*yt);
e=yt-xx*beta;
xxe=xx.*(e*ones(1,length(xx(1,:))));
sehet=sqrt(diag(xxi*xxe'*xxe*xxi));
sehomo=sqrt(diag(xxi*(e'*e))/(ty-n-length(xx(1,:))));
fprintf(out,'Thresholds');
for j=1:length(rrr)
    fprintf(out,'%f   ',rrr(j));
end;
fprintf(out,'\n');
fprintf(out,'\n');
fprintf(out,'Regime-independent Coefficients   standard errors  het standard errors\n');
for j=1:k
    fprintf(out,'%f      %f      %f\n',beta(j),sehomo(j),sehet(j));
end;
fprintf(out,'\n');
fprintf(out,'Regime- dependent Coefficients     standard errors  het standard errors\n');
for j=k+1:k+nr+1
    fprintf(out,'%f      %f      %f\n',beta(j),sehomo(j),sehet(j));
end;
fprintf(out,'\n');
fprintf(out,'\n');

if rep>0
    xx=[xt,ct];
    if abs(max(rr))>0
        j=1;
        while j<=length(rr(:,1))
            xx=[xx,tr(cf.*(thresh<rr(j)))];
            j=j+1;
        end;
    end;
    yp=xx*(yt'/xx')';
    e=yt-yp;
    sse0=e'*e;
    lrt=(sse0/sse1-1)*ty;
    
    fprintf(out,'LR Test for threshold effect:    %f\n',lrt);
    fprintf(out,'\n');
    fprintf(out,'\n');
    stats=zeros(rep,1);
    
    j=1;
    while j<=rep
        eb=reshape(e',tt-1,n)'; %N by T
        %eorig = eb(1:20,1:13)'; %T by N
        ind=0;
        y=eb(ceil(unifrnd(0,1,n,1)*n),:)';  % T by N
        for i=1:length(y(1,:)) %over N
            for jjj=1:length(y(:,1)) % over T
                if ind==0
                    rrrrr=y(jjj,i);
                    ind=1;
                else
                    rrrrr=[rrrrr;y(jjj,i)];
                end;
            end;
        end;
        yb=yp+rrrrr;
        clear rrrrr;
        sse0=sse_calc(yb,[xt,ct]);
        
        [sse1,rhat_b]=r_est(yb,0,trim_);        
        rrr=rhat_b;
        if abs(max(r))>0
            jj=1;
            while jj<=length(r(:,1))
                sse0=sse1;
                [sse1,rhat_b]=r_est(yb,rrr,trim_);
                rrr=[rrr;rhat_b];
                jj=jj+1;
            end;
        end;
        lrt_b=(sse0/sse1-1)*ty;
        stats(j)=lrt_b;
        fprintf(out,'Bootstrap Replication:     %f   %f\n',j,lrt_b);
        j=j+1;
    end;
    fprintf(out,'\n');
    fprintf(out,'\n');
    stats=sortrows(stats,1);
    crits=stats(ceil([.90;.95;.99]*rep));
    fprintf(out,'Number of Bootstrap replications:%f\n',rep);
    fprintf(out,'Bootstrap p-value :              %f\n',mean(stats>lrt));
    fprintf(out,'Critical Values:                 %f\n',crits);
    fprintf(out,'\n');
    fprintf(out,'\n');
    fclose(out);
end;
      
