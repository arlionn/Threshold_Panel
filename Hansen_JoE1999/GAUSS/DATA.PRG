/*   DATA.PRG    */


load invest;
t = 15;
nt = rows(invest);
n = nt/t;

i = invest[.,1];              @ investment/assets                   @
q = invest[.,2];              @ Tobin's Q                           @
c = invest[.,3];              @ cash-flow/assets                    @
d = invest[.,4];              @ debt/assets                         @


output file=data.out reset;

format 12,6;
max_lag = 1;
rhat1 = 0.0157;
rhat2 = 0.5362;

thresh = d;
tt = t-max_lag;

i0 = lag_value(i,n,t,0,max_lag);
q1 = lag_value(q,n,t,1,max_lag);
c1 = lag_value(c,n,t,1,max_lag);
d1 = lag_value(d,n,t,1,max_lag);

xx = i0~q1~c1~d1;
nnt = rows(xx);

j=1; do while j<=4;
  xx[.,j] = sortc(xx[.,j],1);
j=j+1;endo;

qn1=round(nnt/4);
qn2=round(nnt/2);
qn3=round(nnt*.75);
x0 = xx[1,.]';
x1 = xx[qn1,.]';
x2 = xx[qn2,.]';
x3 = xx[qn3,.]';
x4 = xx[nnt,.]';

"Full Sample Summary Statistics";
x0~x1~x2~x3~x4;
"";

e1 = (d1 .<= rhat1);
e2 = (d1 .<= rhat2) - e1;
e3 = 1 - e1 - e2;
f1 = reshape(e1,n,tt);
f2 = reshape(e2,n,tt);
f3 = reshape(e3,n,tt);
g1 = meanc(f1);
g2 = meanc(f2);
g3 = meanc(f3);
g = g1~g2~g3;
g = round(g*100);

st = seqa(1974,1,tt);
"Thresholds " rhat1~rhat2;
"";
"Percentage of Firms in Three Regimes, By Year";
st~g;
"";"";
output off;



proc tr(y,t);
 local yfm,yf;
 yf = reshape(y,rows(y)/t,t);
 yfm = yf - meanc(yf');
 yfm = yfm[.,1:t-1];
retp(vec(yfm'));
endp;


proc lag_value(x,n,t,lagn,max_lag);
local y;
  y = reshape(x,n,t);
  y = y[.,1+max_lag-lagn:t-lagn];
retp(vec(y'));
endp;

