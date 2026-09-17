#include <Rcpp.h>
using namespace Rcpp;

// Scaled recursions with explicit, externally supplied sequence boundaries.
// [[Rcpp::export]]
List fb_sequences_cpp(NumericMatrix emission, NumericVector initial,
                      NumericMatrix transition, LogicalVector start) {
  int n = emission.nrow(), k = emission.ncol();
  NumericMatrix alpha(n,k), beta(n,k), posterior(n,k), counts(k,k);
  NumericVector scale(n), initial_counts(k);
  double ll = 0;
  for (int t=0; t<n; ++t) {
    for (int j=0; j<k; ++j) {
      double pred = start[t] ? initial[j] : 0;
      if (!start[t]) for (int i=0; i<k; ++i) pred += alpha(t-1,i)*transition(i,j);
      alpha(t,j) = pred*emission(t,j); scale[t] += alpha(t,j);
    }
    if (!(scale[t]>0) || !R_finite(scale[t])) stop("Invalid forward scale");
    ll += log(scale[t]);
    for (int j=0; j<k; ++j) alpha(t,j) /= scale[t];
  }
  for (int t=n-1; t>=0; --t) {
    bool last = t==n-1 || start[t+1];
    for (int i=0; i<k; ++i) {
      if (last) beta(t,i)=1;
      else for (int j=0; j<k; ++j)
        beta(t,i) += transition(i,j)*emission(t+1,j)*beta(t+1,j)/scale[t+1];
    }
    double total=0;
    for (int i=0; i<k; ++i) {posterior(t,i)=alpha(t,i)*beta(t,i);total+=posterior(t,i);}
    for (int i=0; i<k; ++i) {
      posterior(t,i)/=total;
      if(start[t]) initial_counts[i]+=posterior(t,i);
    }
    if (!last) {
      double z=0;
      for(int i=0;i<k;++i) for(int j=0;j<k;++j)
        z += alpha(t,i)*transition(i,j)*emission(t+1,j)*beta(t+1,j);
      for(int i=0;i<k;++i) for(int j=0;j<k;++j)
        counts(i,j) += alpha(t,i)*transition(i,j)*emission(t+1,j)*beta(t+1,j)/z;
    }
  }
  return List::create(_["log_likelihood"]=ll, _["posterior"]=posterior,
                      _["filtered"]=alpha, _["transition_sum"]=counts,
                      _["initial_sum"]=initial_counts, _["scales"]=scale);
}

// [[Rcpp::export]]
IntegerVector viterbi_sequences_cpp(NumericMatrix log_emission,
                                    NumericVector initial, NumericMatrix transition,
                                    LogicalVector start) {
  int n=log_emission.nrow(),k=log_emission.ncol();
  NumericMatrix delta(n,k);IntegerMatrix back(n,k);IntegerVector path(n);
  for(int t=0;t<n;++t) for(int j=0;j<k;++j) {
    double best = start[t] ? log(initial[j]) : R_NegInf;
    if (!start[t]) for(int i=0;i<k;++i) {
      double val=delta(t-1,i)+log(transition(i,j));
      if(val>best){best=val;back(t,j)=i;}
    }
    delta(t,j)=best+log_emission(t,j);
  }
  for(int t=n-1;t>=0;--t){
    if(t==n-1 || start[t+1]) {
      double best=R_NegInf;int arg=0;
      for(int j=0;j<k;++j) if(delta(t,j)>best){best=delta(t,j);arg=j;}
      path[t]=arg;
    }else path[t]=back(t+1,path[t+1]);
  }
  for(int t=0;t<n;++t)path[t]+=1;
  return path;
}
