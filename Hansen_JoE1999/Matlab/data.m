% DATA.M %
function data;
load invest.txt;
t=15;
nt=length(invest(:,1));
n=nt/t;

i=invest(:,1);  % investment/assets
q=invest(:,2);  % Tobin's Q 
c=invest(:,3);  % cash-flow/assets
d=invest(:,4);  % debt/assets

% format 12,6; %
max_lag=1;
rhat1=0.0157;
rhat2=0.5362;

thresh=d;
tt=t-max_lag;

i0 = lag_value(i,n,t,0,max_lag);
q1 = lag_value(q,n,t,1,max_lag);
c1 = lag_value(c,n,t,1,max_lag);
d1 = lag_value(d,n,t,1,max_lag);

xx=[i0,q1,c1,d1];

nnt=length(xx(:,1));

j=1;
while j<=4
    xx(:,j)=sortrows(xx(:,j),1);
    j=j+1;
end;

qn1=round(nnt/4);
qn2=round(nnt/2);
qn3=round(nnt*.75);
x0=xx(1,:)';
x1=xx(qn1,:)';
x2=xx(qn2,:)';
x3=xx(qn3,:)';
x4=xx(nnt,:)';

disp('Full Sample Summary Statistics');
for i=1:length(x0);
    fprintf('%f   %f   %f   %f   %f\n',x0(i),x1(i),x2(i),x3(i),x4(i));
end;
disp(' ');

e1=(d1<=rhat1);
e2=(d1<=rhat2)-e1;
e3=1-e1-e2;
f1 = reshape(e1',tt,n)';
f2 = reshape(e2',tt,n)';
f3 = reshape(e3',tt,n)';
g1 = mean(f1)';
g2 = mean(f2)';
g3 = mean(f3)';
g=[g1,g2,g3];
g=round(g*100);

st=1974:(1973+tt);
fprintf('Thresholds:    %f    %f\n',rhat1,rhat2);
disp(' ');
disp('Percentage of Firms in Three Regimes, By Year');
for i=1:length(st)
    fprintf('%f   %f   %f   %f\n',st(i),g(i,1),g(i,2),g(i,3));
end;
disp(' ');
disp(' ');

function r=tr(y,t)
yf=reshape(y',t,length(y(:,q))/t)';
yfm=yf-mean(yf')';
yfm=yfm(:,1:t-1)';
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

function r=lag_value(x,n,t,lagn,max_lag)
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