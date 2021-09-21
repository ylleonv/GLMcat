#include "cdf.h"
#include "reference.h"
#include "adjacentR.h"
#include "sequentialR.h"
#include "cumulativeR.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

// [[Rcpp::export(.GLMcat)]]
List GLMcat(Formula formula,
            DataFrame data,
            std::string ratio,
            Rcpp::List cdf,
            CharacterVector parallel,
            CharacterVector categories_order,
            CharacterVector ref_category,
            std::string threshold,
            Rcpp::List control,
            double normalization){

  // If not cdf given, assume logistic as default
  CharacterVector cdf_given = cdf[0];
  std::string cdf_1;
  if(cdf_given[0] == "NaN"){
    cdf_1 = "logistic";
  }else if(cdf_given[0] != "NaN"){
    std::string cdf2 = cdf[0];
    cdf_1 = cdf2;
  }

  double freedom_degrees = 1;
  double mu = 0;
  if(cdf.size() == 2){
    freedom_degrees = cdf[1];
  }
  if(cdf.size() == 3){
    freedom_degrees = cdf[1];
    mu = cdf[2];
  }

  LogicalVector is_na_ref1 = is_na(categories_order);
  if(is_na_ref1[0]){ // categories order is not given
    categories_order = ref_category; // take ref_cat
  }

  if((ratio != "cumulative") && (threshold != "standard")){
    Rcpp::stop("Unrecognized threshold restriction; for reference, adjacent and sequential ratio the only valid option is standard");
  }

  if (!(threshold == "standard" || threshold == "equidistant" || threshold == "symmetric" )){
    Rcpp::stop("Unrecognized threshold restriction; options are: standard and equidistant");
  }

  class cdf dist1;
  const int N = data.nrows() ; // Number of observations
  // Rcout << parallel << std::endl;

  List Full_M = dist1.All_pre_data_or(formula,
                                      data,
                                      categories_order,
                                      parallel,
                                      threshold,
                                      ratio);

  MatrixXd Y_init = Full_M["Response_EXT"];
  MatrixXd X_EXT = Full_M["Design_Matrix"];
  CharacterVector levs1 = Full_M["Levels"];
  categories_order = Full_M["categories_order"];
  CharacterVector parallel_effect = Full_M["parallel_effect"];
  CharacterVector explanatory_complete = Full_M["Complete_effects"];
  // Rcout << explanatory_complete << std::endl;
  int N_cats = Full_M["N_cats"];

  int P_c = explanatory_complete.length();
  int P_p = 0;
  if(parallel_effect[0] != "NA"){P_p = parallel_effect.length();}

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros or by user option
  Eigen::MatrixXd BETA2;
  BETA2 = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
  Eigen::VectorXd BETA3 = BETA2;

  Eigen::MatrixXd BETA = BETA2;

  NumericVector beta_init;
  double epsilon;
  int iterations_us;

  if(control.size() > 1){

    iterations_us = control[0];
    epsilon = control[1] ;
    beta_init = control[2];

  }else{
    iterations_us = 25;
    epsilon = 1e-06;
    beta_init = R_NaN;
    List control1 = List::create(Named("maxit") = 25 , Named("epsilon") = 1e-06, Named("beta_init")= beta_init);
    control= control1;
    // beta_init = control[2];
  }

  // Rcout << iterations_us << std::endl;
  // cout << epsilon << std::endl;


  // Initialization cumulative special ascending intercepts, it generated some problems
  // if((ratio == "cumulative" && explanatory_complete[0] == "(Intercept)") && threshold == "standard"){
  //   int qm = Q/2;
  //   // Rcout << qm << std::endl;
  //   // qm_vec = std::vector<int> vec(Q-1, qm);
  //   IntegerVector seqvec = seq(1,Q) ; // kind of symmetric around 0
  //   // Rcout << seqvec  << std::endl;
  //   NumericVector seqvec2 = as<NumericVector>(seqvec);
  //   seqvec2 = seqvec2 - qm;
  //   Eigen::Map<Eigen::VectorXd> seqvec1 = as<Eigen::Map<Eigen::VectorXd> >(seqvec2);
  //   // Rcout << seqvec1 << std::endl;
  //   BETA.block(0, 0 , Q , 1) = seqvec1;
  //   Rcout << BETA << std::endl;
  // }

  if(beta_init.length() >= 2 ){
    BETA = as<Eigen::Map<Eigen::VectorXd> >(beta_init);
  }


  std::string Convergence = "False";

  int iteration = 0;
  double Stop_criteria = 1.0;
  MatrixXd X_M_i ;
  VectorXd Y_M_i ;
  VectorXd eta ;
  VectorXd pi ;
  MatrixXd D ;
  MatrixXd Cov_i ;
  MatrixXd W_in ;
  MatrixXd Score_i_2 ;
  MatrixXd F_i_2 ;
  VectorXd LogLikIter;
  LogLikIter = MatrixXd::Zero(1,1) ;
  MatrixXd cov_beta;
  VectorXd Std_Error;
  double LogLik;
  MatrixXd pi_ma(N, Q);
  MatrixXd D_ma(N, Q);
  MatrixXd F_i_final = MatrixXd::Zero(BETA.rows(), BETA.rows());

  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){

  double qp , s0 = 1;

  Environment base_base7("package:base");
  Function my_solve = base_base7["solve"];

  // while ((Stop_criteria>(epsilon / N))){
  while ((Stop_criteria>(epsilon / N)) & (iteration <= (iterations_us ))){
    // Rcout << Stop_criteria << std::endl;

    MatrixXd Score_i = MatrixXd::Zero(BETA.rows(),1);
    MatrixXd F_i = MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;


    // Loop by subject
    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;
      // VectorXd eta1 = eta;
      if(ratio == "reference"){
        ReferenceF ref;
        if(cdf_1 == "logistic"){
          pi = ref.inverse_logistic(eta);
          D = ref.inverse_derivative_logistic(eta);
        }else if(cdf_1 == "normal"){
          pi = ref.inverse_normal(eta);
          D = ref.inverse_derivative_normal(eta);
        }else if(cdf_1 == "cauchy"){
          pi = ref.inverse_cauchy(eta);
          D = ref.inverse_derivative_cauchy(eta);
        }else if(cdf_1 == "gompertz"){
          pi = ref.inverse_gompertz(eta);
          D = ref.inverse_derivative_gompertz(eta);
        }else if(cdf_1 == "gumbel"){
          pi = ref.inverse_gumbel(eta);
          D = ref.inverse_derivative_gumbel(eta);
        }else if(cdf_1 == "laplace"){
          pi = ref.inverse_laplace(eta);
          D = ref.inverse_derivative_laplace(eta);
        }else if(cdf_1 == "student"){
          pi = ref.inverse_student(eta, freedom_degrees);
          D = ref.inverse_derivative_student(eta, freedom_degrees);
        }else if(cdf_1 == "noncentralt"){
          pi = ref.inverse_noncentralt(eta, freedom_degrees, mu);
          D = ref.inverse_derivative_noncentralt(eta, freedom_degrees, mu);
        }else{
          Rcpp::stop("Unrecognized cdf_1; options are: logistic, normal, cauchy, gumbel, gompertz, laplace, student(df), and noncentral(df,mu)");
        }
      }else if(ratio == "adjacent"){
        AdjacentR adj;
        if(cdf_1 == "logistic"){
          pi = adj.inverse_logistic(eta);
          D = adj.inverse_derivative_logistic(eta);
        }else if(cdf_1 == "normal"){
          pi = adj.inverse_normal(eta);
          D = adj.inverse_derivative_normal(eta);
        }else if(cdf_1 == "cauchy"){
          pi = adj.inverse_cauchy(eta);
          D = adj.inverse_derivative_cauchy(eta);
        }else if(cdf_1 == "gompertz"){
          pi = adj.inverse_gompertz(eta);
          D = adj.inverse_derivative_gompertz(eta);
        }else if(cdf_1 == "gumbel"){
          pi = adj.inverse_gumbel(eta);
          D = adj.inverse_derivative_gumbel(eta);
        }else if(cdf_1 == "laplace"){
          pi = adj.inverse_laplace(eta);
          D = adj.inverse_derivative_laplace(eta);
        }else if(cdf_1 == "student"){
          pi = adj.inverse_student(eta, freedom_degrees);
          D = adj.inverse_derivative_student(eta, freedom_degrees);
        }else if(cdf_1 == "noncentralt"){
          pi = adj.inverse_noncentralt(eta, freedom_degrees, mu);
          D = adj.inverse_derivative_noncentralt(eta, freedom_degrees, mu);
        }else{
          Rcpp::stop("Unrecognized cdf_1; options are: logistic, normal, cauchy, gumbel, gompertz, laplace, student(df), and noncentral(df,mu)");
        }
      }else if(ratio == "sequential"){
        SequentialR seq;
        if(cdf_1 == "logistic"){
          pi = seq.inverse_logistic(eta);
          D = seq.inverse_derivative_logistic(eta);
        }else if(cdf_1 == "normal"){
          pi = seq.inverse_normal(eta);
          D = seq.inverse_derivative_normal(eta);
        }else if(cdf_1 == "cauchy"){
          pi = seq.inverse_cauchy(eta);
          D = seq.inverse_derivative_cauchy(eta);
        }else if(cdf_1 == "gompertz"){
          pi = seq.inverse_gompertz(eta);
          D = seq.inverse_derivative_gompertz(eta);
        }else if(cdf_1 == "gumbel"){
          pi = seq.inverse_gumbel(eta);
          D = seq.inverse_derivative_gumbel(eta);
        }else if(cdf_1 == "laplace"){
          pi = seq.inverse_laplace(eta);
          D = seq.inverse_derivative_laplace(eta);
        }else if(cdf_1 == "student"){
          pi = seq.inverse_student(eta, freedom_degrees);
          D = seq.inverse_derivative_student(eta, freedom_degrees);
        }else if(cdf_1 == "noncentralt"){
          pi = seq.inverse_noncentralt(eta, freedom_degrees, mu);
          D = seq.inverse_derivative_noncentralt(eta, freedom_degrees, mu);
        }else{
          Rcpp::stop("Unrecognized cdf_1; options are: logistic, normal, cauchy, gumbel, gompertz, laplace, student(df), and noncentral(df,mu)");
        }
      }else if(ratio == "cumulative"){
        CumulativeR cum;
        if(cdf_1 == "logistic"){
          pi = cum.inverse_logistic(eta);
          D = cum.inverse_derivative_logistic(eta);
        }else if(cdf_1 == "normal"){
          pi = cum.inverse_normal(eta);
          D = cum.inverse_derivative_normal(eta);
        }else if(cdf_1 == "cauchy"){
          pi = cum.inverse_cauchy(eta);
          D = cum.inverse_derivative_cauchy(eta);
        }else if(cdf_1 == "gompertz"){
          pi = cum.inverse_gompertz(eta);
          D = cum.inverse_derivative_gompertz(eta);
        }else if(cdf_1 == "student"){
          pi = cum.inverse_student(eta,freedom_degrees);
          D = cum.inverse_derivative_student(eta,freedom_degrees);
        }else if(cdf_1 == "laplace"){
          pi = cum.inverse_laplace(eta);
          D = cum.inverse_derivative_laplace(eta);
        }else if(cdf_1 == "gumbel"){
          pi = cum.inverse_gumbel(eta);
          D = cum.inverse_derivative_gumbel(eta);
        }else if(cdf_1 == "noncentralt"){
          pi = cum.inverse_noncentralt(eta, freedom_degrees, mu);
          D = cum.inverse_derivative_noncentralt(eta, freedom_degrees, mu);
        }else{
          Rcpp::stop("Unrecognized cdf_1; options are: logistic, normal, cauchy, gumbel, gompertz, laplace, student(df), and noncentral(df,mu)");
        }
      }


      Cov_i = MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      // Rcout << "Cov_i.determinant()" << std::endl;
      // Rcout << Cov_i.determinant() << std::endl;


      // FullPivLU<MatrixXd> lu(Cov_i);
      // bool invertible = lu.isInvertible();
      //
      // if(!invertible || LogLikIter[iteration]>-0.000001 ) {
      //   Rcpp::stop("Fisher matrix is not invertible COV. Check for convergence problems");
      // }

      // NumericMatrix Cov_i_2 = wrap(Cov_i);
      // NumericMatrix Cov_i_3 = my_solve(Cov_i_2, Named("tol")= 2.95222e-200);
      //
      // Eigen::Map<Eigen::MatrixXd> Cov_i_3_4 = as<Eigen::Map<Eigen::MatrixXd> >(Cov_i_3);


      W_in = D * Cov_i.inverse();
      // W_in = D * Cov_i_3_4;
      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );

      // MatrixXd pi_mat = pi;
      // MatrixXd pi_mat1 = pi_mat.transpose();
      // // VectorXd pi_vec1 = pi_mat1;
      //
      // VectorXd pi_vec1(Map<VectorXd>(pi_mat1.data(), pi_mat1.cols()*pi_mat1.rows()));

      pi_ma.row(i) = pi.transpose();
      D_ma.row(i) = D.transpose();
      // pi_ma.row(i) = pi_vec1;
      // pi_ma.row(i) = pi;

      // Rcout << "W_in" << std::endl;
      // Rcout << pi.rows() << std::endl;
      // Rcout << pi_ma.cols() << std::endl;
      // Rcout << pi << std::endl;
      // Rcout << pi_ma << std::endl;

    }

    // VectorXd pima3 = pi_ma.rowwise().sum();
    // pi_ma.resize(pi_ma.rows(), K);
    // VectorXd Ones1 = VectorXd::Ones(pi_ma.rows());
    // pi_ma.col(Q) = Ones1 - pima3 ;

    // To stop when LogLik is smaller than the previous
    if(iteration>10){
      if (LogLikIter[iteration] > LogLik )
        break;
      // iteration = 25;
    }

    //
    // Rcout << iteration << std::endl;
    // Rcout << LogLik << std::endl;

    LogLikIter.conservativeResize(iteration+2, 1);
    LogLikIter(iteration+1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    VectorXd beta_old = BETA;

    // MatrixXd inverse;
    FullPivLU<MatrixXd> lu(F_i);
    bool invertible = lu.isInvertible();

    if(!invertible || LogLik > -0.000001 ) {
      Rcpp::warning("Fisher matrix is not invertible. Check for convergence problems");
    }

    // if(!invertible ) {
    //   Rcpp::stop("Fisher matrix is not invertible. Check for convergence problems");
    // }


    // NumericMatrix F_i_2 = wrap(F_i);
    // NumericMatrix F_i_3 = my_solve(F_i_2, Named("tol")= 2.95222e-200);
    //
    // Eigen::Map<Eigen::MatrixXd> A_eigen = as<Eigen::Map<Eigen::MatrixXd> >(F_i_3);

    BETA = BETA + (F_i.inverse() * Score_i);
    // BETA = BETA + (A_eigen * Score_i);
    // check_tutz = ((BETA - beta_old).norm())/(beta_old.norm()+check_tutz);
    iteration = iteration + 1;
    // if (iteration == 30) {
    //   cout << "Max iter" << endl;
    //   Rcpp::stop("Max iter");
    // }

    F_i_final = F_i;
    // Rcout << "BETA" << std::endl;
    // Rcout << BETA << std::endl;
    // Rcout << "LogLik" << std::endl;
    // Rcout << LogLik << std::endl;


    if(iteration == iterations_us){
      warning("Maximum number of iterations reached.");
    }

  }


  // Rcout << "pi_ma" << std::endl;
  // Rcout << pi_ma << std::endl;
  VectorXd pima3 = pi_ma.rowwise().sum();
  MatrixXd pimadef = pi_ma;
  pimadef.conservativeResize(pi_ma.rows(), K);
  VectorXd Ones1 = VectorXd::Ones(pi_ma.rows());
  pimadef.col(Q) = Ones1 - pima3 ;
  // pimadef.block(0,Q,pi_ma.rows(),1) = Ones1 - pima3 ;

  // Rcout << "pimadef" << std::endl;
  // Rcout << pimadef << std::endl;

  cov_beta = F_i_final.inverse();
  // MatrixXd Hessian = -F_i_final;
  Std_Error = cov_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  std::vector<std::string> text=as<std::vector<std::string> >(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string> >(levs1);

  // Rcout << explanatory_complete << std::endl;
  // Rcout << parallel << std::endl;

  StringVector names;
  if(threshold == "equidistant" || threshold == "symmetric"){
    // if(P_c>1){ // hay alguna complete
    StringVector names1;
    int ind_name = 0;
    if(threshold == "equidistant"){
      StringVector names2(2 + Q*(P_c-1) + P_p);
      for(int var = 0 ; var < 1 ; var++){
        names2[ind_name] = dist1.concatenate(text[var], level_text[0]);
        names2[ind_name+1] = dist1.concatenate(text[var], "distance");
        ind_name = ind_name + 2;
        names1 = names2;
      }}
    else{ // symmetric
      double cath = (Q+1)-((Q+1)/2);

      StringVector names2(cath + Q*(P_c-1) + P_p);

      names2[0] = "center";
      ind_name = ind_name + 1;

      for(int var = 1 ; var < cath ; var++){
        names2[ind_name] = dist1.concatenate("a", level_text[ind_name]);
        // names1[ind_name] = (text[0]);
        ind_name = ind_name + 1;
      }
      names1 = names2;
    }

    if(P_c > 1){
      for(int var = 1 ; var < explanatory_complete.size() ; var++){
        for(int cat = 0 ; cat < Q ; cat++){
          names1[ind_name] = dist1.concatenate(text[var], level_text[cat]);
          ind_name = ind_name + 1;
        }
      }
    }
    if(P_p > 0){
      for(int var_p = 0 ; var_p < P_p ; var_p++){
        names1[ind_name] = parallel_effect[var_p];
        ind_name = ind_name+1;
      }
    }
    names = names1;
  }else{
    StringVector names1(Q*P_c + P_p);
    if(P_c > 0){
      for(int var = 0 ; var < explanatory_complete.size() ; var++){
        for(int cat = 0 ; cat < Q ; cat++){
          names1[(Q*var) + cat] = dist1.concatenate(text[var], level_text[cat]);
        }
      }
    }
    if(P_p > 0){
      for(int var_p = 0 ; var_p < P_p ; var_p++){
        names1[(Q*P_c) + var_p] = parallel_effect[var_p];
      }
    }
    names = names1;
  }

  // Rcout << names<< std::endl;
  // TO NAMED THE RESULT BETAS
  NumericMatrix coef = wrap(BETA);
  rownames(coef) = names;
  int df = coef.length();
  MatrixXd predict_glmcated = X_EXT * BETA;

  VectorXd Ones2 = VectorXd::Ones(Y_init.rows());
  VectorXd vex1 = (Y_init.rowwise().sum()) ;
  Y_init.conservativeResize( Y_init.rows(), K);
  Y_init.col(Q) = (vex1 - Ones2).array().abs() ;
  MatrixXd residuals = Y_init - pimadef;
  VectorXd pi_ma_vec(Map<VectorXd>(pimadef.data(), pimadef.cols()*pimadef.rows()));
  VectorXd Y_init_vec(Map<VectorXd>(Y_init.data(), Y_init.cols()*Y_init.rows()));
  VectorXd div_arr = Y_init_vec.array() / pi_ma_vec.array();
  VectorXd dev_r(Y_init.rows());
  int el_1 = 0;
  for (int element = 0 ; element < div_arr.size() ;  element++){
    if (div_arr[element] != 0){
      dev_r[el_1] = div_arr[element];
      el_1 = el_1 +1 ;
    }
  }
  ArrayXd dev_log = dev_r.array().log();
  double deviance = dev_log.sum();
  deviance = -2*deviance;
  // bool conv = true;


  // Normalization of parameters

  // double qp , s0;

  if((normalization != 1) & (cdf_1 != "logistic")){

    // BETA_3 = BETA_2;
    class Logistic logistic;
    // class Student stu;

    qp = logistic.qdf_logit(normalization);

    if(cdf_1 == "normal"){
      class Normal norm;
      s0 = qp / (norm.qdf_normal(normalization)-norm.qdf_normal(0.5));
    }else if(cdf_1 == "cauchy"){
      class Cauchy cauchy;
      s0 = qp / (cauchy.qdf_cauchy(normalization)- cauchy.qdf_cauchy(0.5));
    }else if(cdf_1 == "gompertz"){
      class Gompertz gompertz;
      s0 = qp / (gompertz.qdf_gompertz(normalization)-gompertz.qdf_gompertz(0.5));
    }else if(cdf_1 == "gumbel"){
      class Gumbel gumbel;
      s0 = qp / (gumbel.qdf_gumbel(normalization)-gumbel.qdf_gumbel(0.5));
    }else if(cdf_1 == "laplace"){
      class Laplace laplace;
      s0 = qp / (laplace.qdf_laplace(normalization)-laplace.qdf_laplace(0.5));
    }else if(cdf_1 == "student"){
      class Student stu;
      s0 = qp / (stu.qdf_student(normalization, freedom_degrees)- stu.qdf_student(0.5, freedom_degrees));
    }else if(cdf_1 == "noncentralt"){
      class Noncentralt noncentralt;
      s0 = qp / (noncentralt.qdf_non_central_t(normalization, freedom_degrees, mu)- noncentralt.qdf_non_central_t(0.5, freedom_degrees, mu));
    }

    // NumericMatrix BETA_3 = BETA_2 * (s0);
    // output_list_dis["normalized_coefficients"] = BETA_3;
  }



  List cdf_list = List::create(
    Named("cdf") = cdf_1,
    Named("freedom_degrees") = freedom_degrees,
    Named("mu") = mu

  );


  if(iteration < iterations_us){
    Convergence = "True";
  }



  List output_list = List::create(
    Named("Function") = "GLMcat",
    Named("coefficients") = coef,
    Named("stderr") = Std_Error,
    Named("iteration") = iteration,
    Named("ratio") = ratio,
    // Named("convergence") = Convergence,
    // Named("exp_variables") = names,
    Named("ref_category") = ref_category,
    Named("threshold") = threshold,
    Named("control") = control,
    Named("normalization_s0") = s0,
    // Named("pinv") = pinv,
    Named("cov_beta") = cov_beta,
    Rcpp::Named("df of the model") = df,
    // Rcpp::Named("fitted") = pi_ma,
    // Rcpp::Named("D_ma") = D_ma,
    // Rcpp::Named("pi_ma_vec") = pi_ma_vec,
    // Rcpp::Named("Y_init_vec") = Y_init_vec,
    // Rcpp::Named("dev_log") = dev_log,
    // Rcpp::Named("deviance") = deviance,
    // Rcpp::Named("residuals") = residuals,
    Named("LogLikelihood") = LogLik,
    // Named("mu_noncentralt") = mu,
    Named("Hessian") = F_i_final,
    Named("LogLikIter") = LogLikIter,
    Named("formula") = formula,
    Named("categories_order") = categories_order,
    Named("parallel") = parallel,
    Named("N_cats") = N_cats,
    Named("nobs_glmcat") = N,
    Named("cdf") = cdf_list
    // Named("coverged") = true,
    // Named("freedom_degrees") = freedom_degrees,
  );

  output_list.attr("class") = "glmcat";

  return output_list;
}




// [[Rcpp::export(.predict_glmcat)]]
NumericMatrix predict_glmcat(List model_object,
                             DataFrame data,
                             String type
){

  std::string ratio = model_object["ratio"];;
  List cdf_list = model_object["cdf"];
  String function = model_object["Function"];
  int N_cats = model_object["N_cats"];
  int N = data.rows();
  CharacterVector names_col;
  List newdataList;
  class cdf dist1;
  Eigen::MatrixXd coef = model_object["coefficients"];

  if(function == "DiscreteCM"){
    ratio = "reference";
    String predict = "predict";
    N = data.rows() / N_cats;
    List arguments = model_object["arguments"];
    names_col = arguments["categories_order"];
    newdataList = dist1.select_data_nested(arguments["formula"],
                                           arguments["case_id"],
                                                    arguments["alternatives"],
                                                             arguments["reference"],
                                                                      arguments["alternative_specific"],
                                                                               data,arguments["intercept"],
                                                                                             predict
    );
  }else{
    names_col = model_object["categories_order"];
    newdataList = dist1.All_pre_data_or(model_object["formula"],
                                        data,
                                        model_object["categories_order"],
                                                    model_object["parallel"],
                                                                model_object["threshold"],
                                                                            model_object["ratio"]);
  }

  MatrixXd X_EXT = newdataList["Design_Matrix"];
  Eigen::MatrixXd predict_glmcat = X_EXT * coef;
  NumericMatrix predict_glmcat1 = wrap(predict_glmcat);
  NumericMatrix predict_glmcat2(N_cats-1, N ,
                                predict_glmcat1.begin());
  NumericMatrix predict_mat = transpose(predict_glmcat2);
  std::string cdf ;

  CharacterVector cdf_given = cdf_list[0];


  std::string cdf_1;
  if(cdf_given[0] == "NaN"){
    cdf = "logistic";
  }else if(cdf_given[0] != "NaN"){
    std::string cdf2 = cdf_list[0];
    cdf = cdf2;
  }
  double freedom_degrees = 1;
  double mu = 0;
  if(cdf_list.size() == 2){
    freedom_degrees = cdf_list[1];
  }
  if(cdf_list.size() == 3){
    freedom_degrees = cdf_list[1];
    mu = cdf_list[2];
  }

  Eigen::VectorXd pi;
  Eigen::MatrixXd pi_total = Eigen::MatrixXd::Zero(N,N_cats-1);
  Eigen::MatrixXd predict_glmcated_eta;
  Eigen::MatrixXd X_M_i;

  for (int i=0; i < N; i++){
    X_M_i = X_EXT.block(i*(N_cats-1) , 0 , N_cats-1 , X_EXT.cols());
    predict_glmcated_eta = X_M_i * coef;
    if(ratio == "reference"){
      ReferenceF ref;
      // Rcout << cdf << std::endl;
      // // print(cdf);
      if(cdf == "logistic"){
        pi = ref.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = ref.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchy"){
        pi = ref.inverse_cauchy(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = ref.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = ref.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = ref.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else if(cdf == "laplace"){
        pi = ref.inverse_laplace(predict_glmcated_eta);
      }else if(cdf == "noncentralt"){
        pi = ref.inverse_noncentralt(predict_glmcated_eta, freedom_degrees, mu);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchy, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "adjacent"){
      AdjacentR adj;
      if(cdf == "logistic"){
        pi = adj.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = adj.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchy"){
        pi = adj.inverse_cauchy(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = adj.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = adj.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = adj.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else if(cdf == "laplace"){
        pi = adj.inverse_laplace(predict_glmcated_eta);
      }else if(cdf == "noncentralt"){
        pi = adj.inverse_noncentralt(predict_glmcated_eta, freedom_degrees, mu);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchy, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "sequential"){
      SequentialR seq;
      // Vector pi depends on selected cdf
      if(cdf == "logistic"){
        pi = seq.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = seq.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchy"){
        pi = seq.inverse_cauchy(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = seq.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = seq.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = seq.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else if(cdf == "laplace"){
        pi = seq.inverse_laplace(predict_glmcated_eta);
      }else if(cdf == "noncentralt"){
        pi = seq.inverse_noncentralt(predict_glmcated_eta, freedom_degrees, mu);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchy, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "cumulative"){
      CumulativeR cum;
      // Vector pi depends on selected cdf
      if(cdf == "logistic"){
        pi = cum.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = cum.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchy"){
        pi = cum.inverse_cauchy(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = cum.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = cum.inverse_student(predict_glmcated_eta,freedom_degrees);
      }else if(cdf == "gumbel"){
        pi = cum.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "laplace"){
        pi = cum.inverse_laplace(predict_glmcated_eta);
      }else if(cdf == "noncentralt"){
        pi = cum.inverse_noncentralt(predict_glmcated_eta, freedom_degrees, mu);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchy, gumbel, gompertz, laplace and student(df)");
      }
    }else{
      Rcpp::stop("Unrecognized radio; options are: reference, adjacent, cumulative and sequential");
    }
    pi_total.row(i) = pi;
  }

  // NumericVector cum_prob = my_rowSums(pi_total);
  NumericVector cum_prob = wrap(pi_total.rowwise().sum());
  // Rcout << cum_prob << std::endl;
  Eigen::Map<Eigen::VectorXd> cum_prob1 = as<Eigen::Map<Eigen::VectorXd> >(cum_prob);
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(pi_total.rows());

  pi_total.conservativeResize(pi_total.rows() , N_cats);
  pi_total.col(N_cats-1) = Ones1 - cum_prob1;

  if(type == "prob"){
    predict_glmcat = pi_total;
    predict_mat = wrap(predict_glmcat);
    colnames(predict_mat) = names_col;
  }else if(type == "linear.predictor"){
    names_col.erase(N_cats-1);
    colnames(predict_mat) = names_col;
  }else{
    Rcpp::stop("Unrecognized type parameter; options are: prob, linear.predictor");
  }
  return predict_mat;
}


// RCPP_MODULE(GLMcatmodule){
//   Rcpp::function("GLMcat", &GLMcat,
//                  List::create(_["formula"],
//                               _["data"],
//                                _["ratio"] = "reference",
//                                _["cdf"] = R_NaN,
//                                _["parallel"] = CharacterVector::create(NA_STRING),
//                                _["categories_order"] = CharacterVector::create(NA_STRING),
//                                _["ref_category"] = CharacterVector::create(NA_STRING),
//                                _["threshold"] = "standard",
//                                _["control"] = R_NaN,
//                                _["normalization"] = 1.0
//                  ),
//                  "GLMcat models");
//
//   Rcpp::function("predict_glmcat", &predict_glmcat,
//                  List::create(_["model_object"] = R_NaN,
//                               _["data"],
//                                _["type"] = "prob"
//                  ),
//                  "GLMcat model predictions");
// }
