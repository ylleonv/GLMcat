#include "cdf.h"
#include "reference.h"
#include "adjacentR.h"
#include "sequentialR.h"
#include "cumulativeR.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

//' Family of models for categorical responses (reference, adjacent, sequential and cumulative ratio)
//'
//' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model.
//' @param ratio a string indicating the F cdf, options are: reference, adjacent, cumulative and sequential. Default value is reference.
//' @param cdf a string indicating the F cdf, options are: logistic, normal, cauchit, student (any df), gompertz, gumbel and laplace.
//' @param categories_order a character vector indicating the incremental order of the categories: c("a", "b", "c"); a<b<c. Alphabetical order is assumed by default. Order is relevant for adjacent, cumulative and sequential ratio.
//' @param ref_category a string indicating the reference category. Proper option for models with reference ratio.
//' @param proportional a character vector indicating the name of the variables with a proportional effect. If variable is categorical, specify the name and the level of the variable as a string "namelevel".
//' @param data a dataframe object in R, with the dependent variable as factor.
//' @param freedom_degrees an optional scalar to indicate the degrees of freedom for the Student cdf.
//' @param threshold restriction to impose on the thresholds, options are: standard, equidistant or symmetric (Valid only for the cumulative ratio).
//' @export
//' @examples
//' data(DisturbedDreams)
//' ref_log_com <- GLMcat(formula = Level ~ Age, data = DisturbedDreams,
//'     ref_category = "Very.severe",
//'     cdf = "logistic", ratio = "reference")
//'
// [[Rcpp::export("GLMcat")]]
List GLMcat(Formula formula,
            DataFrame data,
            std::string ratio,
            std::string cdf,
            CharacterVector proportional,
            CharacterVector categories_order,
            CharacterVector ref_category,
            double freedom_degrees,
            std::string threshold,
            List control){
  // Eigen::VectorXd beta_init){

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
  // Rcout << proportional << std::endl;

  List Full_M = dist1.All_pre_data_or(formula,
                                      data,
                                      categories_order,
                                      proportional,
                                      threshold,
                                      ratio);

  MatrixXd Y_init = Full_M["Response_EXT"];
  MatrixXd X_EXT = Full_M["Design_Matrix"];
  CharacterVector levs1 = Full_M["Levels"];
  categories_order = Full_M["categories_order"];
  CharacterVector explanatory_complete = Full_M["Complete_effects"];
  int N_cats = Full_M["N_cats"];

  int P_c = explanatory_complete.length();
  int P_p = 0;
  if(proportional[0] != "NA"){P_p = proportional.length();}

  int Q = Y_init.cols();
  int K = Q + 1;
  // // // Beta initialization with zeros or by user option
  Eigen::MatrixXd BETA2;
  BETA2 = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
  Eigen::VectorXd BETA3 = BETA2;

  Eigen::MatrixXd BETA = BETA2;
  NumericVector beta_init = control["beta_init"];

  if(beta_init.length() >= 2 ){
    // BETA = beta_init;
    BETA = as<Eigen::Map<Eigen::VectorXd> >(beta_init);
  }
  //
  int iteration = 0;
  // double check_tutz = 1.0;
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
  double epsilon = control["epsilon"] ;
  int iterations_us = control["maxit"];
  while ((Stop_criteria>(epsilon / N)) & (iteration < (iterations_us ))){

    MatrixXd Score_i = MatrixXd::Zero(BETA.rows(),1);
    MatrixXd F_i = MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Loop by subject
    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;
      VectorXd eta1 = eta;
      if(ratio == "reference"){
        ReferenceF ref;
        if(cdf == "logistic"){
          pi = ref.inverse_logistic(eta1);
          D = ref.inverse_derivative_logistic(eta1);
        }else if(cdf == "normal"){
          pi = ref.inverse_normal(eta);
          D = ref.inverse_derivative_normal(eta);
        }else if(cdf == "cauchit"){
          pi = ref.inverse_cauchit(eta);
          D = ref.inverse_derivative_cauchit(eta);
        }else if(cdf == "gompertz"){
          pi = ref.inverse_gompertz(eta);
          D = ref.inverse_derivative_gompertz(eta);
        }else if(cdf == "gumbel"){
          pi = ref.inverse_gumbel(eta);
          D = ref.inverse_derivative_gumbel(eta);
        }else if(cdf == "laplace"){
          pi = ref.inverse_laplace(eta);
          D = ref.inverse_derivative_laplace(eta);
        }else if(cdf == "student"){
          pi = ref.inverse_student(eta, freedom_degrees);
          D = ref.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
        }
      }else if(ratio == "adjacent"){
        AdjacentR adj;
        if(cdf == "logistic"){
          pi = adj.inverse_logistic(eta);
          D = adj.inverse_derivative_logistic(eta);
        }else if(cdf == "normal"){
          pi = adj.inverse_normal(eta);
          D = adj.inverse_derivative_normal(eta);
        }else if(cdf == "cauchit"){
          pi = adj.inverse_cauchit(eta);
          D = adj.inverse_derivative_cauchit(eta);
        }else if(cdf == "gompertz"){
          pi = adj.inverse_gompertz(eta);
          D = adj.inverse_derivative_gompertz(eta);
        }else if(cdf == "gumbel"){
          pi = adj.inverse_gumbel(eta);
          D = adj.inverse_derivative_gumbel(eta);
        }else if(cdf == "laplace"){
          pi = adj.inverse_laplace(eta);
          D = adj.inverse_derivative_laplace(eta);
        }else if(cdf == "student"){
          pi = adj.inverse_student(eta, freedom_degrees);
          D = adj.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
        }
      }else if(ratio == "sequential"){
        SequentialR seq;
        if(cdf == "logistic"){
          pi = seq.inverse_logistic(eta);
          D = seq.inverse_derivative_logistic(eta);
        }else if(cdf == "normal"){
          pi = seq.inverse_normal(eta);
          D = seq.inverse_derivative_normal(eta);
        }else if(cdf == "cauchit"){
          pi = seq.inverse_cauchit(eta);
          D = seq.inverse_derivative_cauchit(eta);
        }else if(cdf == "gompertz"){
          pi = seq.inverse_gompertz(eta);
          D = seq.inverse_derivative_gompertz(eta);
        }else if(cdf == "gumbel"){
          pi = seq.inverse_gumbel(eta);
          D = seq.inverse_derivative_gumbel(eta);
        }else if(cdf == "laplace"){
          pi = seq.inverse_laplace(eta);
          D = seq.inverse_derivative_laplace(eta);
        }else if(cdf == "student"){
          pi = seq.inverse_student(eta, freedom_degrees);
          D = seq.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
        }
      }else if(ratio == "cumulative"){
        CumulativeR cum;
        if(cdf == "logistic"){
          pi = cum.inverse_logistic(eta);
          D = cum.inverse_derivative_logistic(eta);
        }else if(cdf == "normal"){
          pi = cum.inverse_normal(eta);
          D = cum.inverse_derivative_normal(eta);
        }else if(cdf == "cauchit"){
          pi = cum.inverse_cauchit(eta);
          D = cum.inverse_derivative_cauchit(eta);
        }else if(cdf == "gompertz"){
          pi = cum.inverse_gompertz(eta);
          D = cum.inverse_derivative_gompertz(eta);
        }else if(cdf == "student"){
          pi = cum.inverse_student(eta,freedom_degrees);
          D = cum.inverse_derivative_student(eta,freedom_degrees);
        }else if(cdf == "laplace"){
          pi = cum.inverse_laplace(eta);
          D = cum.inverse_derivative_laplace(eta);
        }else if(cdf == "gumbel"){
          pi = cum.inverse_gumbel(eta);
          D = cum.inverse_derivative_gumbel(eta);
        }else{
          Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
        }
      }else{
        Rcpp::stop("Unrecognized radio; options are: reference, adjacent, cumulative and sequential");
      }

      Cov_i = MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      // Rcout << "Cov_i.determinant()" << std::endl;
      // Rcout << Cov_i.determinant() << std::endl;
      W_in = D * Cov_i.inverse();
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
    if(iteration>1){
      if (LogLikIter[iteration] > LogLik)
        break;
      // iteration = 25;
    }
    LogLikIter.conservativeResize(iteration+2, 1);
    LogLikIter(iteration+1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    VectorXd beta_old = BETA;

    // MatrixXd inverse;
    FullPivLU<MatrixXd> lu(F_i);
    bool invertible = lu.isInvertible();

    if(!invertible) {
      Rcpp::stop("Fisher matrix is not invertible");
    }

    BETA = BETA + (F_i.inverse() * Score_i);
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
  Std_Error = cov_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  std::vector<std::string> text=as<std::vector<std::string> >(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string> >(levs1);


  // Rcout << explanatory_complete << std::endl;
  // Rcout << proportional << std::endl;

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

    // print(names1);
    if(P_c > 1){
      for(int var = 1 ; var < explanatory_complete.size() ; var++){
        for(int cat = 0 ; cat < Q ; cat++){
          names1[ind_name] = dist1.concatenate(text[var], level_text[cat]);
          ind_name = ind_name + 1;
        }
      }
    }
    // print(names1);
    if(P_p > 0){
      for(int var_p = 0 ; var_p < proportional.size() ; var_p++){
        names1[ind_name] = proportional[var_p];
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
      for(int var_p = 0 ; var_p < proportional.size() ; var_p++){
        names1[(Q*P_c) + var_p] = proportional[var_p];
      }
    }
    names = names1;
  }

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

  List output_list = List::create(
    Named("coefficients") = coef,
    Named("stderr") = Std_Error,
    Named("iteration") = iteration,
    Named("ratio") = ratio,
    // Named("AIC") = AIC,
    // Named("pinv") = pinv,
    Named("cov_beta") = cov_beta,
    Rcpp::Named("df of the model") = df,
    Rcpp::Named("fitted") = pi_ma,
    Rcpp::Named("D_ma") = D_ma,
    // Rcpp::Named("pi_ma_vec") = pi_ma_vec,
    // Rcpp::Named("Y_init_vec") = Y_init_vec,
    // Rcpp::Named("dev_log") = dev_log,
    Rcpp::Named("deviance") = deviance,
    // Rcpp::Named("residuals") = residuals,
    Named("LogLikelihood") = LogLik,
    // Named("freedom_degrees") = freedom_degrees,
    // Named("Y_init") = Y_init,
    // Named("LogLikIter") = LogLikIter,
    Named("formula") = formula,
    Named("categories_order") = categories_order,
    Named("proportional") = proportional,
    Named("N_cats") = N_cats,
    Named("nobs_glmcat") = N,
    Named("cdf") = cdf,
    Named("coverged") = true,
    Named("freedom_degrees") = freedom_degrees
  );

  output_list.attr("class") = "glmcat";

  return output_list;
}



//' GLMcat model predictions
//'
//' @param model_object a GLMcat model
//' @param data a data frame with the predictor variables used in the GLMcat model.
//' @param type The type of prediction to obtain. \code{"prob"} gives probabilities,
//' \code{"cum.prob"} gives cumulative probabilities and \code{"linear.predict"} gives
//' the linear predictor.
//' @rdname predict_glmcat
//' @export
//' @examples
//' data(DisturbedDreams)
//' mod1 <- GLMcat(formula = Level ~ Age,
//' data = DisturbedDreams, cdf = "logistic")
//' predict_glmcat(mod1, data = DisturbedDreams[1:5, ], type = "prob")
// [[Rcpp::export("predict_glmcat")]]
NumericMatrix predict_glmcat(List model_object,
                             DataFrame data,
                             String type
){

  class cdf dist1;
  // Environment base_env("package:base");
  // Function my_rowSums = base_env["rowSums"];
  int N_cats = model_object["N_cats"];
  Eigen::MatrixXd coef = model_object["coefficients"];

  // model_object["coefficients"]

  List newdataList = dist1.All_pre_data_NEWDATA(model_object["formula"],
                                                data,
                                                model_object["categories_order"],
                                                            model_object["proportional"],
                                                                        N_cats);

  Eigen::MatrixXd Design_Matrix = newdataList["Design_Matrix"];
  Eigen::MatrixXd predict_glmcated_eta;

  std::string ratio = model_object["ratio"];
  String cdf = model_object["cdf"];
  double freedom_degrees = model_object["freedom_degrees"];

  Eigen::VectorXd pi;
  int N = data.rows();
  Eigen::MatrixXd X_M_i;
  Eigen::MatrixXd pi_total = Eigen::MatrixXd::Zero(N,N_cats-1);
  for (int i=0; i < N; i++){

    X_M_i = Design_Matrix.block(i*(N_cats-1) , 0 , N_cats-1 , Design_Matrix.cols());
    predict_glmcated_eta = X_M_i * coef;

    if(ratio == "reference"){
      ReferenceF ref;
      if(cdf == "logistic"){
        pi = ref.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = ref.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchit"){
        pi = ref.inverse_cauchit(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = ref.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = ref.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = ref.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "adjacent"){
      AdjacentR adj;
      if(cdf == "logistic"){
        pi = adj.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = adj.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchit"){
        pi = adj.inverse_cauchit(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = adj.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = adj.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = adj.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "sequential"){
      SequentialR seq;
      // Vector pi depends on selected cdf
      if(cdf == "logistic"){
        pi = seq.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = seq.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchit"){
        pi = seq.inverse_cauchit(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = seq.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "gumbel"){
        pi = seq.inverse_gumbel(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = seq.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
      }
    }else if(ratio == "cumulative"){
      CumulativeR cum;
      // Vector pi depends on selected cdf
      if(cdf == "logistic"){
        pi = cum.inverse_logistic(predict_glmcated_eta);
      }else if(cdf == "normal"){
        pi = cum.inverse_normal(predict_glmcated_eta);
      }else if(cdf == "cauchit"){
        pi = cum.inverse_cauchit(predict_glmcated_eta);
      }else if(cdf == "gompertz"){
        pi = cum.inverse_gompertz(predict_glmcated_eta);
      }else if(cdf == "student"){
        pi = cum.inverse_student(predict_glmcated_eta,freedom_degrees);
      }else if(cdf == "gumbel"){
        pi = cum.inverse_gumbel(predict_glmcated_eta);
      }else{
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchit, gumbel, gompertz, laplace and student(df)");
      }
    }else{
      Rcpp::stop("Unrecognized radio; options are: reference, adjacent, cumulative and sequential");
    }
    pi_total.row(i) = pi;
  }

  // NumericVector cum_prob = my_rowSums(pi_total);
  NumericVector cum_prob = wrap(pi_total.rowwise().sum());
  Eigen::Map<Eigen::VectorXd> cum_prob1 = as<Eigen::Map<Eigen::VectorXd> >(cum_prob);
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(pi_total.rows());

  pi_total.conservativeResize(pi_total.rows() , N_cats);
  pi_total.col(N_cats-1) = Ones1 - cum_prob1;

  NumericVector predict_glmcat;
  NumericMatrix predict_mat;
  CharacterVector names_col = model_object["categories_order"];

  if(type == "prob"){
    predict_glmcat = pi_total;
    predict_mat = wrap(predict_glmcat);
    colnames(predict_mat) = names_col;
  }else if(type == "linear.predict"){
    names_col.erase(N_cats-1);
    predict_glmcat = Design_Matrix * coef;
    NumericMatrix predict_glmcat1( N_cats-1, pi_total.rows() , predict_glmcat.begin());
    predict_mat = transpose(predict_glmcat1);
    colnames(predict_mat) = names_col;
  }else{
    Rcpp::stop("Unrecognized type parameter; options are: prob, linear.predict");
  }
  return predict_mat;
}


RCPP_MODULE(GLMcatmodule){
  Rcpp::function("GLMcat", &GLMcat,
                 List::create(_["formula"],
                              _["data"],
                               _["ratio"] = "reference",
                               _["cdf"] = "logistic",
                               _["proportional"] = CharacterVector::create(NA_STRING),
                               _["categories_order"] = CharacterVector::create(NA_STRING),
                               _["ref_category"] = CharacterVector::create(NA_STRING),
                               _["freedom_degrees"] = R_NaN,
                               _["threshold"] = "standard",
                               _["control"] = List::create(_["maxit"] = 25, _["epsilon"] = 1e-07, _["beta_init"] = NumericVector::create(NA_REAL))
                               // _["beta_init"] = NumericVector::create(NA_REAL)
                                // _["beta_init"] = R_NaN
                 ),
                 "GLMcat models");
  Rcpp::function("predict_glmcat", &predict_glmcat,
                 List::create(_["model_object"] = R_NaN,
                              _["data"],
                               _["type"] = "prob"
                 ),
                 "GLMcat model predictions");
}
