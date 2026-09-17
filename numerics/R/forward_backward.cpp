#include <Rcpp.h>
using namespace Rcpp;

// Scaled recursions identical to the reference R implementation.
// [[Rcpp::export]]
List forward_backward_compiled(NumericMatrix emission, NumericVector initial,
                               NumericMatrix transition) {
  int T = emission.nrow(), K = emission.ncol();
  NumericMatrix alpha(T,K), beta(T,K), posterior(T,K), counts(K,K);
  NumericVector scale(T);
  double ll = 0;
  for (int t=0; t<T; ++t) {
    for (int j=0; j<K; ++j) {
      double pred = t == 0 ? initial[j] : 0;
      if (t > 0) for (int i=0; i<K; ++i) pred += alpha(t-1,i)*transition(i,j);
      alpha(t,j)=pred*emission(t,j); scale[t]+=alpha(t,j);
    }
    if (!(scale[t]>0) || !R_finite(scale[t])) stop("Invalid forward scale");
    ll += std::log(scale[t]);
    for (int j=0; j<K; ++j) alpha(t,j)/=scale[t];
  }
  for (int j=0; j<K; ++j) beta(T-1,j)=1;
  for (int t=T-2; t>=0; --t) for (int i=0; i<K; ++i) {
    for (int j=0; j<K; ++j) beta(t,i)+=transition(i,j)*emission(t+1,j)*beta(t+1,j);
    beta(t,i)/=scale[t+1];
  }
  for (int t=0; t<T; ++t) {
    double total=0;
    for (int j=0; j<K; ++j) {posterior(t,j)=alpha(t,j)*beta(t,j);total+=posterior(t,j);}
    for (int j=0; j<K; ++j) posterior(t,j)/=total;
    if (t<T-1) {
      double denom=0;
      for (int i=0; i<K; ++i) for (int j=0; j<K; ++j)
        denom+=alpha(t,i)*transition(i,j)*emission(t+1,j)*beta(t+1,j);
      for (int i=0; i<K; ++i) for (int j=0; j<K; ++j)
        counts(i,j)+=alpha(t,i)*transition(i,j)*emission(t+1,j)*beta(t+1,j)/denom;
    }
  }
  return List::create(_["log_likelihood"]=ll, _["posterior"]=posterior,
                      _["transition_sum"]=counts);
}
