/* 
This is the C++ file to be compiled and loaded by Rcpp (in R/RStudio). 
To use it, typically no modification is required. 
*/

/* # include <Rcpp.h>
 */
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* here we are given that x are distinct and have already been sorted,
 * while y corresponds to x are also given 
 * return the fitted values for all xs as a vector
 */
List fLmTwoCasts(NumericMatrix Xr, NumericVector yr) {
  int n = Xr.nrow(), k = Xr.ncol();
  arma::mat X(Xr.begin(), n, k, false);
  arma::colvec y(yr.begin(), yr.size(), false);
  // int df = n - k;
  
  // fit model y ~ X, extract residuals
  arma::colvec coef = arma::solve(X, y);
  arma::colvec res  = y - X*coef;
  arma::colvec fitted  = X*coef;
  return List::create(coef, res, fitted);
}



/*
 * inner product between two vectors
 */
double innerprod(NumericVector x, NumericVector y) {
  int n = x.size();
  double sum = 0;
  for (int i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}

/*
 * inner product between a vector and a group of basis functions of shape _/,
 * the basis are given b, y correspond to their sorted x values
 */
NumericVector innerprod_seq(NumericVector y, NumericVector b){
  int n = y.size();
  NumericVector res(n);
  res[n-1]=0;
  double nor = 1/sqrt(min(diff(b)));
  double sum_y_right = y[n-1];
  for (int i = n-2; i >=0 ; i--){
    res[i]=res[i+1]+nor*(b[i+1]-b[i])*sum_y_right;
    sum_y_right += y[i];
  }
  return res/nor;
}

/*
 * given information from fitting previous x_1,...,x_n-1
 * compute new results for fitting x_1,...x_n
 * List contains three elements, first is the knot location from previous iter
 * (not necessarily sorted), 
 * NOTE - the first element is always zero -- imaginary knot location for the constant function
 * second is the corresponding coefficients for basis functions, with the first 
 * being that for constant
 * third is the fitted values for f(x_1),...,f(x_n) (so must be sorted)
 */
List nextinfo(List previnfo, NumericVector x, NumericVector y){
  NumericVector knot = previnfo[0];
  NumericVector coeff = previnfo[1];
  NumericVector fitted = previnfo[2];
  int n = y.size();
  int i, j, k;
    
  fitted.push_back(y[n-1]);
  
  double y_thres = (sum(coeff)-coeff[0])*(x[n-1]-x[n-2])+fitted[n-2];
  
  if(y[n-1] >= y_thres){
    if(y[n-1] > y_thres){
      knot.push_back(x[n-2]);
      coeff.push_back((y[n-1]-y_thres)/(x[n-1]-x[n-2]));
    }
    return List::create(knot, coeff, fitted);
  }
  else {
 
    // first, creat a vector of (0,...,0,1)
    NumericVector yind(n);
    yind[n-1]=1;
    // now create a vector of y but with the last element = the threshold
    NumericVector yatthres = clone(y);
    fitted[n-1] =  y_thres;
    
    while(true){
      
      yatthres[n-1]=y_thres;
      
      // then create required design matrix
      NumericMatrix dm (n, knot.size());
      for (k = 0; k < n; k++) dm(k,0)=1;
      for (j = 1; j < knot.size(); j++)
        for (k = 0; k < n; k++)
          dm(k,j) = (x[k] - knot[j] > 0) ? (x[k] - knot[j]) : 0;
      // then run linear regression on that
      List res1 = fLmTwoCasts(dm, yind);
      NumericVector coeffyind= res1[0];
      NumericVector residyind= res1[1];
      NumericVector fittedyind= res1[2];
      NumericVector temp = coeff/coeffyind;
      //check the prime feasibility

      double primeStep = 0;
      int rm_index = -1;

      if (temp.size()>1){
        for(i = temp.size()-1; i > 0; i--){
          if (coeffyind[i] > 1e-6 && temp[i] > 1e-9){
            if( (primeStep > temp[i]) || (rm_index == -1)){
              rm_index = i;
              primeStep = temp[rm_index];
            }
          }
        }
      }
      // check the dual feasibility
      List res2 = fLmTwoCasts(dm, yatthres);
      NumericVector part1 = innerprod_seq(res2[1],x);
      NumericVector part2 = -innerprod_seq(res1[1],x);
      
      for (j=0; j < part2.size(); j++) if (fabs(part2[j])< 1e-8) part2[j] = 1e-8 * (part2[j]>0?1:-1);
      NumericVector temp2 = -part1/part2;
  
      double dualStep = 0;
      int add_index = -1;
      if (temp2.size()>=1){
        // note that temp2 has first index representing the knot at the smallest observation (i.e. linear function)
        // while knot has that represent the constant knot
        // so please be careful about their indices...
        for(i = 0; i <  temp2.size()-1; i++){
          if (temp2[i] > 1e-6 && (part1[i]) < -1e-6 && (part2[i]) > 1e-6){
          // also check that index to be added cannot be in (or too close to) the existing knot set
            if (knot.size()>=2 && min(abs(x[i]-knot[Range(1,knot.size()-1)])) > 1e-6 ){
              if((temp2[i] < dualStep) || (add_index == -1)){
                add_index = i;
                dualStep = temp2[i];
              }
            }
          }
        }
      }
      
      //decide add or subtract
      // nothing to remove, so (can be shown) that must be constant
      if (rm_index ==-1){
        coeff = coeff- (y_thres-y[n-1]) * coeffyind;
        fitted = fitted - (y_thres-y[n-1])  * fittedyind;
        return List::create(knot, coeff, fitted);
      }
      else {
        if((primeStep < dualStep) || add_index == -1){
          // job done - 
          if (y_thres - primeStep <= y[n-1]){
            coeff = coeff- (y_thres-y[n-1]) * coeffyind;
            fitted = fitted - (y_thres-y[n-1])  * fittedyind;
            return List::create(knot, coeff, fitted);
          } 
          // not yet done, subtract index
          else{
            coeff = coeff- primeStep * coeffyind;
            fitted = fitted - primeStep * fittedyind;
            knot.erase(rm_index);
            coeff.erase(rm_index);
            y_thres = y_thres - primeStep;
          }
        }
        else if ((primeStep > dualStep)){
          // job done - 
          if (y_thres - dualStep <= y[n-1]){
            coeff = coeff- (y_thres-y[n-1]) * coeffyind;            
            fitted = fitted - (y_thres-y[n-1])  * fittedyind;
            return List::create(knot, coeff, fitted);
          } 
          // not yet done, add index
          else{
            coeff = coeff- dualStep * coeffyind;
            fitted = fitted - dualStep * fittedyind;
            knot.push_back(x[add_index]);
            coeff.push_back(0.0);
            y_thres = y_thres - dualStep;
          }
        }
        else{
          // must be primeStep == dualStep
          // job done - 
          if (y_thres - primeStep <= y[n-1]){
            coeff = coeff- (y_thres-y[n-1]) * coeffyind;
            fitted = fitted - (y_thres-y[n-1])  * fittedyind;
            return List::create(knot, coeff, fitted);
          } 
          // subtract and add
          else{
            coeff = coeff- primeStep * coeffyind;
            fitted = fitted - primeStep * fittedyind;
            knot.erase(rm_index);
            coeff.erase(rm_index);
            knot.push_back(x[add_index]);
            coeff.push_back(0.0);
            y_thres = y_thres - primeStep;
          }
        }
      }
    }
  }
}


/* Update the coefficients and the fitted values in previnfo
* to avoid numerical inaccuracy of double type from accumulating.
* So far not used [not really necessary after controlling the 
* epsilons in the main code to be 1e-6]
*/
List restable(List thisinfo, NumericVector x, NumericVector y){
  NumericVector knot = thisinfo[0];
  int n = y.size();
  int j, k;
  // then create required design matrix
  NumericMatrix dm (n, knot.size());
  for (k = 0; k < n; k++) dm(k,0)=1;
  for (j = 1; j < knot.size(); j++)
    for (k = 0; k < n; k++)
      dm(k,j) = (x[k] - knot[j] > 0) ? (x[k] - knot[j]) : 0;
  // then run linear regression on that
  List res1 = fLmTwoCasts(dm, y);
  NumericVector coeff= res1[0];
  NumericVector fitted= res1[2];
  return List::create(knot, coeff, fitted);
}  

/* Order the elements of y by sorting x
*/
NumericVector Rcpp_sort(NumericVector x, NumericVector y) {
  // First create a vector of indices
  IntegerVector idy = seq_along(y) - 1;
  // Then sort that vector by the values of y
  std::sort(idy.begin(), idy.end(), [&](int i, int j){return x[i] < x[j];});
  // And return x in that order
  return y[idy];
}

/* reverse the order of elements in the list 
*/
NumericVector rcppRev(NumericVector x) {
    NumericVector revX = clone<NumericVector>(x);
    std::reverse(revX.begin(), revX.end());
    ::Rf_copyMostAttrib(x, revX); 
    return revX;
}
  


/* cvxreg implements (increasing and) convex sequential algorithm 
* (with the increasing monotonicity constraint)
*/
// [[Rcpp::export]]
List cvxreg(NumericVector x, NumericVector y) {
  int n = x.size();
  int i;
  
  // sort the observations according to x, if not yet done so
  y = Rcpp_sort(x,y);
  x = Rcpp_sort(x,x);
  
  //initialization with one observation
  NumericVector knot(1);
  NumericVector coeff(1);
  NumericVector fitted(1);
  
  knot[0]=0; coeff[0]=y[0]; fitted[0]=y[0];
  List previnfo = List::create(knot, coeff, fitted);
  for(i = 1; i < n; i++){
    NumericVector xx = x[Range(0,i)];
    NumericVector yy = y[Range(0,i)];
    previnfo = nextinfo(previnfo, xx, yy);
    // restable might be needed once a while? not in this version
    // previnfo = restable(previnfo, xx, yy);
  }
  return previnfo;
}
  
/* sshapedreg implements S-shaped regression using cvxreg 
* (with the increasing monotonicity constraint)
*/  
// [[Rcpp::export]]
List sshapedreg(NumericVector x, NumericVector y) {
  int n = x.size();
  int i, j, k, inflex_index;
  double rss, inflex;
  // vector with a single constant knot, to be only used for initialization here
  // for knot_all_l and knot_all_r
  NumericVector knot_constant(1);
  
  // sort the observations according to x, if not yet done so
  y = Rcpp_sort(x,y);
  x = Rcpp_sort(x,x);
  
  NumericVector xrev = -rcppRev(x); 
  NumericVector yrev = -rcppRev(y);
  
  
  NumericVector fitted(n);
  
  NumericVector rss_l(n);
  NumericVector f_l_inflection(n);
  NumericVector df_l_inflection(n);
  List knot_all_l = List::create(knot_constant);
  NumericVector rss_r(n);
  NumericVector f_r_inflection(n);
  NumericVector df_r_inflection(n);  
  List knot_all_r = List::create(knot_constant);
  
  // first pass left to right convex increasing
  //initialization with one observation
  NumericVector knot_l(1),knot_r(1);
  NumericVector coeff_l(1),coeff_r(1);
  NumericVector fitted_l(1),fitted_r(1);
  
  knot_l[0]=0; coeff_l[0]=y[0]; fitted_l[0]=y[0]; f_l_inflection[0]=y[0];
  List previnfo_l = List::create(knot_l, coeff_l, fitted_l);
  for(i = 1; i < n; i++){
    NumericVector xx_l = x[Range(0,i)];
    NumericVector yy_l = y[Range(0,i)];
    previnfo_l = nextinfo(previnfo_l, xx_l, yy_l);
    NumericVector newknot_l = previnfo_l[0];
    NumericVector newfitted_l = previnfo_l[2];
    rss_l[i] = sum((newfitted_l-yy_l)*(newfitted_l-yy_l));
    f_l_inflection[i]=newfitted_l[i];
    df_l_inflection[i]=sum(coeff_l) - coeff_l[0];
    knot_all_l.push_back(newknot_l);
  }
  
  // second pass right to left concave increasing, so apply cvxreg algorithm 
  // using -x and -y (i.e xrev and yrev)
  knot_r[0]=0; coeff_r[0]=yrev[0]; fitted_r[0]=yrev[0]; f_r_inflection[0]=yrev[0];
  List previnfo_r = List::create(knot_r, coeff_r, fitted_r);
  for(i = 1; i < n; i++){
    NumericVector xx_r = xrev[Range(0,i)];
    NumericVector yy_r = yrev[Range(0,i)];
    previnfo_r = nextinfo(previnfo_r, xx_r, yy_r);
    NumericVector newknot_r = previnfo_r[0];
    NumericVector newfitted_r = previnfo_r[2];
    rss_r[i] = sum((newfitted_r-yy_r)*(newfitted_r-yy_r));
    f_r_inflection[i]=newfitted_r[i];
    df_r_inflection[i] = sum(coeff_r) - coeff_r[0];
    knot_all_r.push_back(newknot_r);
  }
  
  
  inflex_index = 0; rss = rss_r[n-1]; inflex = x[inflex_index];
  // here in the algorithm inflection points belongs to the right subinterval
  // i.e. [0,...,inflex_index-1] [inflex_index,...,n-1]
  for(i = 1; i < n; i++){
    if(x[i]-x[i-1] > 1e-7){
      double slope_mid = (-f_r_inflection[n-1-i] - f_l_inflection[i-1])/(x[i]-x[i-1]);
      if ((df_l_inflection[i-1] <= slope_mid) && (df_r_inflection[n-1-i] <= slope_mid) && (rss_l[i-1] + rss_r[n-1-i]< rss)){
        inflex_index = i;
        rss = rss_l[i-1] + rss_r[n-1-i];
        inflex = x[inflex_index];
      }
    }
    else if ((f_l_inflection[i-1] <= -f_r_inflection[n-1-i]) && (rss_l[i-1] + rss_r[n-1-i]< rss)){
      inflex_index = i;
      rss = rss_l[i-1] + rss_r[n-1-i];
      inflex = x[inflex_index];
    }
  }
  if (rss_l[n-1] < rss){
    inflex_index = n;
    rss = rss_l[n-1];
    inflex = x[n-1];
  }
  // concave
  if (inflex_index == 0){
    NumericVector knot = knot_all_r[n-1];
    // then create required design matrix
    NumericMatrix dm (n, knot.size());
    for (k = 0; k < n; k++) dm(k,0)=1;
    for (j = 1; j < knot.size(); j++)
      for (k = 0; k < n; k++)
        dm(k,j) = (xrev[k] - knot[j] > 0) ? (xrev[k] - knot[j]) : 0;
    // then run linear regression on that
    List res1 = fLmTwoCasts(dm, yrev);
    fitted = -rcppRev(res1[2]);
  }
  //convex
  else if (inflex_index == n){
    NumericVector knot = knot_all_l[n-1];
    // then create required design matrix
    NumericMatrix dm (n, knot.size());
    for (k = 0; k < n; k++) dm(k,0)=1;
    for (j = 1; j < knot.size(); j++)
      for (k = 0; k < n; k++)
        dm(k,j) = (x[k] - knot[j] > 0) ? (x[k] - knot[j]) : 0;
    // then run linear regression on that
    List res1 = fLmTwoCasts(dm, y);
    fitted = res1[2];
  }
  // proper S-shape
  else{
    NumericVector knot_l = knot_all_l[inflex_index-1];
    NumericVector xx_l = x[Range(0,inflex_index-1)];
    NumericVector yy_l = y[Range(0,inflex_index-1)];
    int n_l = inflex_index;
    
    
    // then create required design matrix
    NumericMatrix dm_l (n_l, knot_l.size());
    for (k = 0; k < n_l; k++) dm_l(k,0)=1;
    for (j = 1; j < knot_l.size(); j++)
      for (k = 0; k < n_l; k++)
        dm_l(k,j) = (xx_l[k] - knot_l[j] > 0) ? (xx_l[k] - knot_l[j]) : 0;
    // then run linear regression on that
    List res1 = fLmTwoCasts(dm_l, yy_l);
    NumericVector fitted_l = res1[2];
    for (i = 0; i < n_l; i++) fitted[i]= fitted_l[i];
    
    NumericVector knot_r = knot_all_r[n-inflex_index-1];
    NumericVector xx_r = xrev[Range(0, n-inflex_index-1)];
    NumericVector yy_r = yrev[Range(0, n-inflex_index-1)];
    int n_r = n - inflex_index;
    
    NumericMatrix dm_r (n_r, knot_r.size());
    for (k = 0; k < n_r; k++) dm_r(k,0)=1;
    for (j = 1; j < knot_r.size(); j++)
      for (k = 0; k < n_r; k++)
        dm_r(k,j) = (xx_r[k] - knot_r[j] > 0) ? (xx_r[k] - knot_r[j]) : 0;
    // then run linear regression on that
    List res2 = fLmTwoCasts(dm_r, yy_r);
    NumericVector fitted_r = -rcppRev(res2[2]);
    for (i = 0; i < n_r; i++) fitted[i+n_l]= fitted_r[i];
  }
  
  return List::create(inflex, fitted, rss);
}
