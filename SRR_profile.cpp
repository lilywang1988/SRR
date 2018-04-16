# include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
List SRR_profile_cpp(vec Y,vec Facility, mat z, vec sizev,double criterion,int max_iter,float bound) {
  int N, nF;
  N=Y.size();
  vec Facility_list=unique(Facility);
  nF=Facility_list.size();
  vec gamma_est(nF);
  gamma_est.fill(0);
  int beta_length = z.n_cols;
  vec beta_est(beta_length);
  beta_est.fill(0);
  int rep_n = 0;
  vec llk_all(max_iter);
  vec beta_update(beta_length);
  beta_update.fill(1000);
  vec Y_p(N);
  vec gamma_update_pt1(nF);
  vec gamma_update_pt2(nF);
  vec gamma_update(nF);
  vec beta_U(z.n_cols);
  mat beta_I(z.n_cols,z.n_cols);
  do{rep_n+=1;
    cout<<rep_n<<endl;
    vec gamma_est_long(N);
    int start=0;
    for(int f=0; f<nF; f++){
      int end=start+sizev(f)-1;
      gamma_est_long(span(start,end))=vec(sizev(f)).fill(gamma_est(f));
      start=end+1;
    }
    vec exp_eta=exp(gamma_est_long+z*beta_est); 
    vec p = exp_eta/(exp_eta+1); // / is element-wise division
    double llk = sum(log(p));
    llk_all(rep_n-1) = llk;
    vec q = 1-p;
    vec p_times_q = p%q;//%element-wise multiplication

    vec gamma_update_pt1(nF);
    gamma_update_pt1.fill(0);

    gamma_update_pt2.fill(0);
    Y_p=Y-p;
    start=0;
    for(int f=0; f<nF; f++){
      int end=start+sizev(f)-1;
      gamma_update_pt1[f]=sum(Y_p(span(start,end)));
      gamma_update_pt2[f]=sum(p_times_q(span(start,end)));
      start=end+1;
    }
    gamma_update=gamma_update_pt1/gamma_update_pt2;
    gamma_est = gamma_est + gamma_update;
    start=0;
    for(int f=0; f<nF; f++){
      int end=start+sizev(f)-1;
      gamma_est_long(span(start,end))=vec(sizev(f)).fill(gamma_est(f));
      start=end+1;
    }
    exp_eta=exp(gamma_est_long+z*beta_est);
    p = exp_eta/(exp_eta+1);
    p_times_q = p%q;
    Y_p=Y-p;
    vec beta_U = z.t()*Y_p;
    mat z_pq(size(z));
    for(unsigned int c=0;c<z.n_cols;c++){
      z_pq.col(c)=z.col(c)%p_times_q;
    }
    mat beta_I = z.t()* z_pq;
    beta_update = inv(beta_I)*beta_U;//inv() is to inverse a matrix
    beta_est = beta_est +  beta_update;
    uvec gamma_index1=find(gamma_est > bound);
    uvec gamma_index2=find(gamma_est < -bound);
    gamma_est.elem(gamma_index1) = vec(gamma_index1.size()).fill(bound);
    gamma_est.elem(gamma_index2) = vec(gamma_index2.size()).fill(-bound);
    uvec beta_index1=find(beta_est > bound);
    uvec beta_index2=find(beta_est < -bound);
    beta_est.elem(beta_index1) = vec(beta_index1.size()).fill(bound);
    beta_est.elem(beta_index2) = vec(beta_index2.size()).fill(-bound);
  } while (rep_n < (max_iter-1) && max(abs((beta_update))) > criterion);
  llk_all.resize(rep_n);
  List res;
  res["beta.est"] = beta_est;
  res["gamma.est"] = gamma_est;
  res["llk_all"] = llk_all;
  res["rep.n"] = rep_n-1;
  return(res);
}
