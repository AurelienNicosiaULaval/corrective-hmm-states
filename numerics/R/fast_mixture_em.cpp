#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Full scaled EM, with the same safeguards as mixture_em_loop in R.
// Sequence starts have a common, freely estimated initial distribution.
// [[Rcpp::export]]
List fast_mixture_em_cpp(NumericVector y, NumericVector init0,
 NumericMatrix trans0, NumericMatrix w0, NumericMatrix mu0, NumericMatrix sd0,
 IntegerVector counts, LogicalVector start, int max_iter=2000,
 double tolerance=1e-8, double sigma_min=.05) {
 int n=y.size(), k=init0.size(), m=w0.ncol();
 if(n<1 || !start[0] || start.size()!=n) stop("Invalid sequence boundaries");
 NumericVector initial=clone(init0); NumericMatrix tr=clone(trans0),w=clone(w0),mu=clone(mu0),sd=clone(sd0);
 std::vector<int> state,component;
 for(int i=0;i<k;i++) for(int j=0;j<counts[i];j++){state.push_back(i);component.push_back(j);}
 int r=state.size(), ns=0;for(int t=0;t<n;t++)ns+=start[t];
 std::vector<double> em(n*k),cd(n*r),a(n*k),b(n*k),post(n*k),scale(n),tc(k*k),ic(k),tot(r),sy(r),sv(r),trace;
 bool finish=false,converged=false; int updates=0;
 for(int iteration=0;iteration<=max_iter;iteration++) {
  if(iteration%50==0)checkUserInterrupt();
  std::fill(em.begin(),em.end(),0);
  for(int h=0;h<r;h++) {
   int i=state[h],j=component[h];
   for(int t=0;t<n;t++) {double d=w(i,j)*std::max(1e-300,R::dnorm(y[t],mu(i,j),sd(i,j),false));cd[t*r+h]=d;em[t*k+i]+=d;}
  }
  double ll=0;
  for(int t=0;t<n;t++) {
   double z=0;
   for(int j=0;j<k;j++) {
    em[t*k+j]=std::max(1e-300,em[t*k+j]);
    double pred=start[t]?initial[j]:0;
    if(!start[t])for(int i=0;i<k;i++)pred+=a[(t-1)*k+i]*tr(i,j);
    a[t*k+j]=pred*em[t*k+j];z+=a[t*k+j];
   }
   if(!(z>0) || !R_finite(z))stop("Invalid forward scale");
   scale[t]=z;ll+=log(z);for(int j=0;j<k;j++)a[t*k+j]/=z;
  }
  std::fill(tc.begin(),tc.end(),0);std::fill(ic.begin(),ic.end(),0);
  for(int t=n-1;t>=0;t--) {
   bool last=t==n-1 || start[t+1];double z=0;
   for(int i=0;i<k;i++) {
    double v=last?1:0;
    if(!last)for(int j=0;j<k;j++)v+=tr(i,j)*em[(t+1)*k+j]*b[(t+1)*k+j]/scale[t+1];
    b[t*k+i]=v;post[t*k+i]=a[t*k+i]*v;z+=post[t*k+i];
   }
   for(int i=0;i<k;i++){post[t*k+i]/=z;if(start[t])ic[i]+=post[t*k+i];}
   if(!last) {
    double zt=0;for(int i=0;i<k;i++)for(int j=0;j<k;j++)zt+=a[t*k+i]*tr(i,j)*em[(t+1)*k+j]*b[(t+1)*k+j];
    for(int i=0;i<k;i++)for(int j=0;j<k;j++)tc[i*k+j]+=a[t*k+i]*tr(i,j)*em[(t+1)*k+j]*b[(t+1)*k+j]/zt;
   }
  }
  if(!trace.empty() && ll<trace.back()-1e-6*(1+std::abs(trace.back())))stop("EM likelihood decreased");
  converged=!trace.empty() && ll-trace.back()<tolerance*(1+std::abs(trace.back()));
  trace.push_back(ll);
  if(finish || iteration==max_iter)break;
  double iz=0;for(int i=0;i<k;i++){initial[i]=ic[i]/ns+1e-10;iz+=initial[i];}
  for(int i=0;i<k;i++){
   initial[i]/=iz;double z=0;for(int j=0;j<k;j++){tr(i,j)=tc[i*k+j]+1e-10;z+=tr(i,j);}for(int j=0;j<k;j++)tr(i,j)/=z;
  }
  std::fill(tot.begin(),tot.end(),0);std::fill(sy.begin(),sy.end(),0);std::fill(sv.begin(),sv.end(),0);
  for(int t=0;t<n;t++)for(int h=0;h<r;h++){
   int i=state[h];double resp=post[t*k+i]*cd[t*r+h]/em[t*k+i];
   cd[t*r+h]=resp;tot[h]+=resp;sy[h]+=resp*y[t];
  }
  for(int h=0;h<r;h++)mu(state[h],component[h])=sy[h]/(tot[h]+1e-12);
  for(int t=0;t<n;t++)for(int h=0;h<r;h++){double d=y[t]-mu(state[h],component[h]);sv[h]+=cd[t*r+h]*d*d;}
  for(int h=0;h<r;h++){int i=state[h],j=component[h];w(i,j)=tot[h]+1e-12;sd(i,j)=sqrt(std::max(sigma_min*sigma_min,sv[h]/(tot[h]+1e-12)));}
  for(int i=0;i<k;i++){
   double z=0;for(int j=0;j<counts[i];j++)z+=w(i,j);for(int j=0;j<counts[i];j++)w(i,j)/=z;
   std::vector<int> order(counts[i]);for(int j=0;j<counts[i];j++)order[j]=j;
   std::sort(order.begin(),order.end(),[&](int u,int v){return mu(i,u)<mu(i,v);});
   std::vector<double> ww(counts[i]),mm(counts[i]),ss(counts[i]);
   for(int j=0;j<counts[i];j++){ww[j]=w(i,order[j]);mm[j]=mu(i,order[j]);ss[j]=sd(i,order[j]);}
   for(int j=0;j<counts[i];j++){w(i,j)=ww[j];mu(i,j)=mm[j];sd(i,j)=ss[j];}
  }
  updates++;if(converged && iteration>0)finish=true;
 }
 NumericMatrix posterior(n,k);for(int t=0;t<n;t++)for(int i=0;i<k;i++)posterior(t,i)=post[t*k+i];
 return List::create(_["initial"]=initial,_["transition"]=tr,_["weights"]=w,_["means"]=mu,_["sds"]=sd,
 _["posterior"]=posterior,_["log_likelihood"]=trace.back(),_["log_likelihood_trace"]=wrap(trace),
 _["n_iter"]=updates,_["converged"]=converged,_["component_counts"]=counts,_["K"]=k,_["M"]=m);
}
