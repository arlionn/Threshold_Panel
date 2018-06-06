/*
THRESH_P.PRG

This is a GAUSS program file.
It replicates the estimation, testing and graphs reported in
"Threshold Effects in Non-Dynamic Panels:
Estimation, Testing and Inference"

For questions, please contact

Bruce E. Hansen
Department of Economics
Social Science Building
University of Wisconsin
Madison, WI 53706-1393
bhansen@ssc.wisc.edu
http://www.ssc.wisc.edu/~bhansen/


This program file loads the GAUSS dataset "invest.fmt".
It creates the output file "thresh.out"

*/
library pgraph;
load invest;
t = 15;
nt = rows(invest);
n = nt/t;

i = invest[.,1];     @ investment/assets                             @
q = invest[.,2];     @ Tobin's Q                                     @
c = invest[.,3];     @ cash-flow/assets                              @
d = invest[.,4];     @ debt/assets                                   @


qn = 400;            @ number of quantiles to examine                @
conf_lev = .95;      @ confidence level for threshold                @
_vgraph = 1;         @ set to 1 to graph likelihood ratios           @
_boot_1 = 3;       @ # of replications, 0 for no bootstrap, single (300) @
_boot_2 = 3;       @ # of replications, 0 for no bootstrap, double (300) @
_boot_3 = 3;       @ # of replications, 0 for no bootstrap, triple (300) @
_trim_1 = .01;       @ percentage to trim before search, single      @
_trim_2 = .01;       @ percentage to trim before search, double      @
_trim_3 = .05;       @ percentage to trim before search, triple      @

output file=thresh.out reset;  output off;

max_lag = 1;
tt = t-max_lag;
ty = n*(t-max_lag-1);

y  = lag_v(i,0);  yt = tr(y);
cf = lag_v(c,1);  ct = tr(cf);
q1 = lag_v(q,1);
d1 = lag_v(d,1);      @ set to threshold variable @

x = q1~(q1.^2)~(q1.^3)~d1~(q1.*d1);
k = cols(x);
xt = zeros(rows(yt),k);
j=1; do while j<=k;
  xt[.,j] = tr(x[.,j]);
j=j+1;endo;
thresh = d1;
dd = unique(thresh,1);
qnt1 = qn*_trim_1;
sq = seqa(_trim_1,1/qn,qn-2*qnt1+1);
qq1 = dd[floor(sq*rows(dd))];
qn1 = rows(qq1);
cc = -2*ln(1-sqrt(conf_lev));

output on;
"Number of Firms        " n;
"Number of Years used   " tt;
"Total Observations     " ty;
"Number of Quantiles    " qn;
"Confidence Level       " conf_lev;
"";"";
"*******************************************************";
"";"";
sse0 = sse_calc(yt,xt~ct);
"Zero Threshold Model";
"Sum of Squared Errors                   " sse0;
"";"";
"*******************************************************";
"";"";

"Single Threshold Model";
"";
output off;
rhat1 = model(0,_trim_1,_boot_1,0);
output on;
"*******************************************************";
"";"";

"Double Threshold Model";
"Trimming Percentage    " _trim_2;
"";
"First Iteration";
output off;
rhat2 = model(rhat1,_trim_2,_boot_2,2);
output on;
"Second Iteration";
output off;
rhat1 = model(rhat2,_trim_2,0,1);
output on;
"";"";
"*******************************************************";
"";"";

"Triple Threshold Model";
"Trimming Percentage    " _trim_3;
"";
output off;
rhat3 = model(rhat1|rhat2,_trim_3,_boot_3,3);
output on;
"";"";
"*******************************************************";
"";"";
output off;


/****************** PROCS ************************/

proc tr(y);
 local yfm,yf;
 yf = reshape(y,n,tt);
 yfm = yf - meanc(yf');
 yfm = yfm[.,1:tt-1];
retp(vec(yfm'));
endp;

proc lag_v(x,lagn);
local y;
  y = reshape(x,n,t);
  y = y[.,1+max_lag-lagn:t-lagn];
retp(vec(y'));
endp;

proc sse_calc(y,x);
local e;
  e = y - x*(y/x);
retp(e'e);
endp;


proc thr_sse(y,q,r);
local n,sse,qi,rr,d,xx;
n = rows(q);
sse = zeros(n,1);
qi = 1; do while qi<=n;
  if r==0; rr = q[qi]; else; rr = r|q[qi]; endif;
  rr = sortc(rr,1);
  xx = xt~ct;
  j = 1; do while j <= rows(rr);
     d = (thresh .< rr[j]);
     xx = xx~tr(cf.*d);
  j=j+1;endo;
  sse[qi] = sse_calc(y,xx);
qi = qi+1;endo;
retp(sse);
endp;


proc (2) = r_est(y,r,_trim);
local qq,rr,i,nn,ii,sse,rihat,qnt;
if maxc(r) .== 0;
  qq = qq1; rr = 0;
else;
  rr = sortc(r,1);
  i = seqa(1,1,qn1);
  nn = sumc(qq1 .< (rr'))';
  qnt = qn*_trim;
  ii = ((i .<= (nn+qnt))-(i .<= (nn-qnt)))*ones(rows(rr),1);
  qq = delif(qq1,ii);
endif;
sse = thr_sse(y,qq,rr);
rihat = minindc(sse);
retp(sse[rihat],qq[rihat]);
endp;


proc (1) = model(r,_trim,rep,it);
local qq,rr,i,nn,ii,sse,rihat,rhat,sse1,lr,rhats,xx,crits,
rrr,e,sse0,yp,lrt,stats,j,eb,yb,yp_b,lrt_b,jj,dd,d,xxi,beta,
sehet,sehomo,titname,tit,xname,xlab,nr,rhat_b,qnt;

if maxc(r)==0;
  qq = qq1; rr = 0;
else;
  rr = sortc(r,1);
  i = seqa(1,1,qn1);
  nn = sumc(qq1 .< (rr'));
  qnt = qn*_trim;
  ii=((i.<=(nn+qnt)')-(i.<=(nn-qnt)'))*ones(rows(rr),1);
  qq = delif(qq1,ii);
endif;
sse = thr_sse(yt,qq,rr);
rihat = minindc(sse);
rhat = qq[rihat];
sse1 = sse[rihat];
lr = (sse/sse1 - 1)*ty;
rhats = selif(qq,(lr .< cc));
if _vgraph==1;
  graphset;
  if     it==0; 
    titname="Figure 1\LConfidence Interval Construction in Single Threshold Model";
    xname="Threshold Parameter";
  elseif it==1; 
    titname = "Figure 3\LConfidence Interval Construction in Double Threshold Model"; 
    xname="First Threshold Parameter";
  elseif it==2; 
    titname = "Figure 2\LConfidence Interval Construction in Double Threshold Model"; 
    xname="Second Threshold Parameter";
  elseif it==3; 
    titname = "Confidence Interval Construction in Triple Threshold Model"; 
    xname="Third Threshold Parameter";
  endif;
  _pdate="";
  ylabel("Likelihood Ratio");
  title(titname);
  xlabel(xname);
  xy(qq,lr~(ones(rows(qq),1)*cc));
endif;
output on;
if maxc(r) .ne 0;
  "Fixed Thresholds       " rr';
  rrr = sortc((rr|rhat),1);
else;
  rrr = rhat;
endif;
"Threshold Estimate     " rhat;
"Confidence Region      " minc(rhats)~maxc(rhats);
"Sum of Squared Errors  " sse1;
"Trimming Percentage    " _trim;
"";"";
nr = rows(rrr);
xx = xt;
dd = zeros(rows(thresh),nr);
j=1; do while j<=nr;
  dd[.,j] = (thresh .< rrr[j]);
  d = dd[.,j];
  if j>1;
    d = d - dd[.,j-1];
  endif;
  xx = xx~tr(cf.*d);
j=j+1;endo;
d = 1-dd[.,nr];
xx = xx~tr(cf.*d);
xxi = invpd(moment(xx,0));
beta = xxi*(xx'yt);
e = yt - xx*beta;
sehet = sqrt(diag(xxi*moment(xx.*e,0)*xxi));
sehomo = sqrt(diag(xxi*(e'e)/(ty-n-cols(xx))));
beta = beta~sehomo~sehet;
"Thresholds";
rrr';
"";
"Regime-independent Coefficients, standard errors, het standard errors";
beta[1:k,.];
"";
"Regime-dependent Coefficients, standard errors, het standard errors";
beta[k+1:k+nr+1,.];
"";"";
output off;
if rep .> 0;
  xx = xt~ct;
  if maxc(rr) .ne 0;
    j = 1; do while j <= rows(rr);
       xx = xx~tr(cf.*(thresh .< rr[j]));
    j=j+1;endo;
  endif;
  yp = xx*(yt/xx);
  e = yt-yp;
  sse0 = e'e;
  lrt = (sse0/sse1-1)*ty;
  output on;
  "LR Test for threshold effect  " lrt;
  output off;
  "";"";
  stats = zeros(rep,1);
  j=1; do while j<=rep;
    eb   = reshape(e,n,tt-1);
    yb   = yp + vec(eb[ceil(rndu(n,1)*n),.]');
    sse0 = sse_calc(yb,xt~ct);
    {sse1,rhat_b} = r_est(yb,0,_trim);
    rrr = rhat_b;
    if maxc(r) .ne 0;
      "rows(r)=";rows(r);
      jj = 1; do while jj<=rows(r);
        sse0 = sse1;
        {sse1,rhat_b} = r_est(yb,rrr,_trim);
        rrr = rrr|rhat_b;
      jj = jj + 1; endo;
      jj;
    endif;
    lrt_b = (sse0/sse1-1)*ty;
    stats[j] = lrt_b;
    "Bootstrap Replication " j~lrt_b;
  j=j+1;endo;
  "";"";
  stats = sortc(stats,1);
  crits = stats[ceil((.90|.95|.99)*rep)];
  output on;
  "Number of Bootstrap replications   " rep;
  "Bootstrap p-value                  " meanc(stats .> lrt);
  "Critical Values   " crits;
  "";"";
  output off;
endif;
retp(rhat);
endp;


/***************************************************/

