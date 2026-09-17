#include <Rcpp.h>
using namespace Rcpp;

// Interval likelihood for rounded movement rates and optional turning angles.
// All transitions are row-stochastic; the last category is the logit reference.
NumericVector probabilities(const NumericVector &p, int &offset, int n) {
  NumericVector z(n); double maxz=0;
  for(int j=0;j<n-1;j++){z[j]=p[offset++];maxz=std::max(maxz,z[j]);}
  double total=0;
  for(int j=0;j<n;j++){z[j]=std::exp(z[j]-maxz);total+=z[j];}
  return z/total;
}
double log_difference(double a, double b) {
  if(!R_finite(b)) return a;
  if(b>=a) return R_NegInf;
  return a+std::log(-std::expm1(b-a));
}
double movement_log_cdf(double x,double loc,double shape,int family,bool lower) {
  if(x<=0) return lower ? R_NegInf : 0.0;
  if(family==0) return R::pgamma(x,std::exp(shape),std::exp(loc-shape),lower,1);
  if(family==1) return R::pweibull(x,std::exp(shape),std::exp(loc),lower,1);
  if(family==2) return R::pnorm(std::log(x),loc,std::exp(shape),lower,1);
  return R::pnorm(std::log1p(x),loc,std::exp(shape),lower,1);
}
double interval_mass(double lo,double hi,double loc,double scale,int family) {
  double upper=movement_log_cdf(hi,loc,scale,family,true);
  double result;
  if(upper>-0.6931471805599453) {
    double a,b;
    // For the truncated normal, an interval beginning at zero starts at Y=0.
    if(family==3 && lo<=0) a=R::pnorm(0,loc,std::exp(scale),false,true);
    else a=movement_log_cdf(lo,loc,scale,family,false);
    b=movement_log_cdf(hi,loc,scale,family,false);
    result=log_difference(a,b);
  } else {
    double lower=(family==3 && lo<=0) ? R::pnorm(0,loc,std::exp(scale),true,true) : movement_log_cdf(lo,loc,scale,family,true);
    result=log_difference(upper,lower);
  }
  if(family==3) result-=R::pnorm(0,loc,std::exp(scale),false,true);
  return std::max(result,-745.0);
}

// [[Rcpp::export]]
List movement_model_cpp(NumericVector par,NumericVector lower,NumericVector upper,
                        NumericVector angle,IntegerVector counts,int family,
                        bool joint=false,bool detail=false,bool shared_angle=false) {
  int K=counts.size(),N=lower.size(),R=0,offset=0;
  for(int k=0;k<K;k++)R+=counts[k];
  NumericVector initial=probabilities(par,offset,K);
  NumericMatrix transition(K,K);
  for(int k=0;k<K;k++){
    NumericVector row=probabilities(par,offset,K);
    for(int j=0;j<K;j++)transition(k,j)=row[j];
  }
  NumericVector weight(R),location(R),scale(R),direction(R),rho(R);
  IntegerVector group(R);int r=0;
  for(int k=0;k<K;k++){
    NumericVector w=probabilities(par,offset,counts[k]);
    for(int a=0;a<counts[k];a++,r++){weight[r]=w[a];group[r]=k;}
  }
  for(int a=0;a<R;a++){location[a]=par[offset++];scale[a]=par[offset++];}
  if(joint){
    if(shared_angle){
      for(int k=0;k<K;k++){
        double mu=par[offset++],concentration=par[offset++];
        for(int a=0;a<R;a++)if(group[a]==k){direction[a]=mu;rho[a]=concentration;}
      }
    }else for(int a=0;a<R;a++){direction[a]=par[offset++];rho[a]=par[offset++];}
  }
  if(offset!=par.size())stop("Parameter length mismatch");
  NumericMatrix log_components(N,R),log_emission(N,K);
  for(int t=0;t<N;t++){
    for(int a=0;a<R;a++){
      double v=R_IsNA(lower[t]) ? 0.0 : interval_mass(lower[t],upper[t],location[a],scale[a],family);
      if(joint && R_finite(angle[t]))v+=std::log1p(-rho[a]*rho[a])-std::log(2*M_PI)-std::log(1+rho[a]*rho[a]-2*rho[a]*std::cos(angle[t]-direction[a]));
      log_components(t,a)=v;
    }
    for(int k=0;k<K;k++){
      double maxv=R_NegInf;
      for(int a=0;a<R;a++)if(group[a]==k)maxv=std::max(maxv,std::log(weight[a])+log_components(t,a));
      double s=0;
      for(int a=0;a<R;a++)if(group[a]==k)s+=std::exp(std::log(weight[a])+log_components(t,a)-maxv);
      log_emission(t,k)=maxv+std::log(s);
    }
  }
  NumericVector predicted=clone(initial),filtered(K),next(K),scores(N);
  double ll=0;
  for(int t=0;t<N;t++){
    double maxv=R_NegInf;
    for(int k=0;k<K;k++)maxv=std::max(maxv,log_emission(t,k));
    double total=0;
    for(int k=0;k<K;k++){filtered[k]=predicted[k]*std::exp(log_emission(t,k)-maxv);total+=filtered[k];}
    scores[t]=maxv+std::log(total);ll+=scores[t];
    for(int j=0;j<K;j++){next[j]=0;for(int k=0;k<K;k++)next[j]+=filtered[k]/total*transition(k,j);}
    predicted=clone(next);
  }
  if(!detail)return List::create(_["nll"]=-ll);
  return List::create(_["nll"]=-ll,_["log_score"]=scores,_["initial"]=initial,
    _["transition"]=transition,_["weight"]=weight,_["location"]=location,
    _["scale"]=scale,_["direction"]=direction,_["rho"]=rho,_["group"]=group+1,
    _["log_emission"]=log_emission,_["log_components"]=log_components);
}

// [[Rcpp::export]]
List movement_water_cpp(NumericVector par, NumericVector lower, NumericVector upper,
                        NumericVector angle, int family, NumericVector water,
                        bool detail=false) {
  int P=par.size(),N=lower.size();
  NumericVector base(P-2);
  for(int j=0;j<P-2;j++)base[j]=par[j];
  List emissions=movement_model_cpp(base,lower,upper,angle,IntegerVector::create(1,1),family,true,true,false);
  NumericMatrix loge=emissions["log_emission"];
  NumericVector predicted=clone(as<NumericVector>(emissions["initial"])),scores(N);
  NumericMatrix filtered(N,2),p_first(N,2);
  double ll=0;
  for(int t=0;t<N;t++) {
    if(t>0) {
      for(int k=0;k<2;k++)p_first(t,k)=R::plogis(base[1+k]+par[P-2+k]*water[t],0,1,true,false);
      predicted[0]=filtered(t-1,0)*p_first(t,0)+filtered(t-1,1)*p_first(t,1);
      predicted[1]=1-predicted[0];
    }
    double maxv=std::max(loge(t,0),loge(t,1)),total=0;
    for(int k=0;k<2;k++){filtered(t,k)=predicted[k]*std::exp(loge(t,k)-maxv);total+=filtered(t,k);}
    scores[t]=maxv+std::log(total);ll+=scores[t];
    for(int k=0;k<2;k++)filtered(t,k)/=total;
  }
  if(!detail)return List::create(_["nll"]=-ll);
  return List::create(_["nll"]=-ll,_["log_score"]=scores,_["filtered"]=filtered,
                      _["p_first"]=p_first,_["emissions"]=emissions);
}
