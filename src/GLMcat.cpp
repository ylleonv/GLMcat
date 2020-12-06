#include "distribution.h"
#include "reference.h"
#include "adjacentR.h"
#include "sequentialR.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

//' Family of models for categorical responses (reference, adjacent and sequential ratio)
//'
//' @param formula a symbolic description of the model to be fit. An expression of the form y ~ model is interpreted as a specification that the response y is modelled by a linear predict_glmcator specified symbolically by model.
//' @param ratio an string indicating the F distribution, options are: reference, adjacent, cumulative and sequential.
//' @param distribution an string indicating the F distribution, options are: logistic, normal, cauchit, student (any df), gompertz, gumbel.
//' @param categories_order a character vector indicating the incremental order of the categories: c("a", "b", "c"); a<b<c. Alphabetical order is assumed by default.
//' @param proportional a character vector indicating the name of the variables with a proportional effect.
//' @param data a dataframe object in R, with the dependent variable as factor.
//' @param freedom_degrees an optional scalar to indicate the degrees of freedom for the Student distribution.
//' @return GLMcat returns a list which can be examined with the function summary.
//' @export
//' @examples
//' data(DisturbedDreams)
//' GLMcat(formula = Level ~ Age, data = DisturbedDreams,
//' distribution = "logistic",ratio = "reference")
// [[Rcpp::export("GLMcat")]]
List GLMcat(Formula formula,
            std::string ratio, std::string distribution,
            CharacterVector categories_order,
            CharacterVector proportional,
            DataFrame data,
            double freedom_degrees){

  // S4 x("pcglm");

  class distribution dist1;
  const int N = data.nrows() ; // Number of observations
  List Full_M = dist1.All_pre_data_or(formula, data,
                                      categories_order, proportional, ratio);

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
  // // // Beta initialization with zeros
  MatrixXd BETA;
  BETA = MatrixXd::Zero(X_EXT.cols(),1);
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
  MatrixXd var_beta;
  VectorXd Std_Error;
  double LogLik;
  MatrixXd pi_ma(N, K);
  MatrixXd F_i_final = MatrixXd::Zero(BETA.rows(), BETA.rows());

  // for (int iteration=1; iteration < 18; iteration++){
  // while (check_tutz > 0.0001){
  double epsilon = 0.0001 ;
  while ((Stop_criteria>(epsilon / N)) & (iteration < 26)){

    MatrixXd Score_i = MatrixXd::Zero(BETA.rows(),1);
    MatrixXd F_i = MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Loop by subject
    for (int i=0; i < N; i++){
      // Block of size (p,q), starting at (i,j): matrix.block(i,j,p,q);
      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      Y_M_i = Y_init.row(i);
      eta = X_M_i * BETA;

      if(ratio == "reference"){
        ReferenceF ref;
        if(distribution == "logistic"){
          pi = ref.inverse_logistic(eta);
          D = ref.inverse_derivative_logistic(eta);
        }else if(distribution == "normal"){
          pi = ref.inverse_normal(eta);
          D = ref.inverse_derivative_normal(eta);
        }else if(distribution == "cauchit"){
          pi = ref.inverse_cauchit(eta);
          D = ref.inverse_derivative_cauchit(eta);
        }else if(distribution == "gompertz"){
          pi = ref.inverse_gompertz(eta);
          D = ref.inverse_derivative_gompertz(eta);
        }else if(distribution == "gumbel"){
          pi = ref.inverse_gumbel(eta);
          D = ref.inverse_derivative_gumbel(eta);
        }else if(distribution == "student"){
          pi = ref.inverse_student(eta, freedom_degrees);
          D = ref.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
        }
      }else if(ratio == "adjacent"){
        AdjacentR adj;
        if(distribution == "logistic"){
          pi = adj.inverse_logistic(eta);
          D = adj.inverse_derivative_logistic(eta);
        }else if(distribution == "normal"){
          pi = adj.inverse_normal(eta);
          D = adj.inverse_derivative_normal(eta);
        }else if(distribution == "cauchit"){
          pi = adj.inverse_cauchit(eta);
          D = adj.inverse_derivative_cauchit(eta);
        }else if(distribution == "gompertz"){
          pi = adj.inverse_gompertz(eta);
          D = adj.inverse_derivative_gompertz(eta);
        }else if(distribution == "gumbel"){
          pi = adj.inverse_gumbel(eta);
          D = adj.inverse_derivative_gumbel(eta);
        }else if(distribution == "student"){
          pi = adj.inverse_student(eta, freedom_degrees);
          D = adj.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
        }
      }else if(ratio == "sequential"){
        SequentialR seq;
        // Vector pi depends on selected distribution
        if(distribution == "logistic"){
          pi = seq.inverse_logistic(eta);
          D = seq.inverse_derivative_logistic(eta);
        }else if(distribution == "normal"){
          pi = seq.inverse_normal(eta);
          D = seq.inverse_derivative_normal(eta);
        }else if(distribution == "cauchit"){
          pi = seq.inverse_cauchit(eta);
          D = seq.inverse_derivative_cauchit(eta);
        }else if(distribution == "gompertz"){
          pi = seq.inverse_gompertz(eta);
          D = seq.inverse_derivative_gompertz(eta);
        }else if(distribution == "gumbel"){
          pi = seq.inverse_gumbel(eta);
          D = seq.inverse_derivative_gumbel(eta);
        }else if(distribution == "student"){
          pi = seq.inverse_student(eta, freedom_degrees);
          D = seq.inverse_derivative_student(eta, freedom_degrees);
        }else{
          Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
        }
      }else{
        Rcpp::stop("Unrecognized radio; options are: reference, adjacent, cumulative and sequential");
      }

      Cov_i = MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());

      // Rcout << "Cov_i.determinant()" << std::endl;
      // Rcout << Cov_i.determinant() << std::endl;
      W_in = D * Cov_i.inverse();
      // Rcout << "W_in" << std::endl;
      // Rcout << W_in << std::endl;


      Score_i_2 = X_M_i.transpose() * W_in * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );

      pi_ma.row(i) = pi.transpose();

    }

    VectorXd Ones1 = VectorXd::Ones(pi_ma.rows());
    pi_ma.col(Q) = Ones1 - pi_ma.rowwise().sum() ;

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

  // var_beta = (((X_EXT.transpose() * F_i_final) * X_EXT).inverse());


  CompleteOrthogonalDecomposition<MatrixXd> cqr(F_i_final);
  MatrixXd pinv = cqr.pseudoInverse();

  var_beta = F_i_final.inverse();
  Std_Error = var_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  std::vector<std::string> text=as<std::vector<std::string>>(explanatory_complete);
  std::vector<std::string> level_text=as<std::vector<std::string>>(levs1);
  StringVector names(Q*P_c + P_p);
  if(P_c > 0){
    for(int var = 0 ; var < explanatory_complete.size() ; var++){
      for(int cat = 0 ; cat < Q ; cat++){
        names[(Q*var) + cat] = dist1.concatenate(text[var], level_text[cat]);
      }
    }
  }
  if(P_p > 0){
    for(int var_p = 0 ; var_p < proportional.size() ; var_p++){
      names[(Q*P_c) + var_p] = proportional[var_p];
    }
  }

  // TO NAMED THE RESULT BETAS
  NumericMatrix coef = wrap(BETA);
  rownames(coef) = names;

  int df = coef.length();

  // IC
  // double AIC = (-2*LogLik) + (2 *coef.length());
  // double BIC = (-2*LogLik) + (coef.length() * log(N) );

  MatrixXd predict_glmcated = X_EXT * BETA;

  VectorXd Ones2 = VectorXd::Ones(Y_init.rows());
  VectorXd vex1 = (Y_init.rowwise().sum()) ;
  Y_init.conservativeResize( Y_init.rows(), K);
  Y_init.col(Q) = (vex1 - Ones2).array().abs() ;
  MatrixXd residuals = Y_init - pi_ma;
  VectorXd pi_ma_vec(Map<VectorXd>(pi_ma.data(), pi_ma.cols()*pi_ma.rows()));
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

  List output_list = List::create(
    Named("coefficients") = coef,
    Named("stderr") = Std_Error,
    Named("iteration") = iteration,
    Named("ratio") = ratio,
    // Named("AIC") = AIC,
    Named("pinv") = pinv,
    Named("var_beta") = var_beta,
    Rcpp::Named("df of the model") = df,
    // Rcpp::Named("predict_glmcated") = predict_glmcated,
    // Rcpp::Named("fitted") = pi_ma,
    // Rcpp::Named("pi_ma_vec") = pi_ma_vec,
    // Rcpp::Named("Y_init_vec") = Y_init_vec,
    // Rcpp::Named("dev_log") = dev_log,
    Rcpp::Named("deviance") = deviance,
    // Rcpp::Named("residuals") = residuals,
    Named("Log-likelihood") = LogLik,
    // Named("freedom_degrees") = freedom_degrees,
    // Named("Y_init") = Y_init,
    // Named("LogLikIter") = LogLikIter,
    Named("formula") = formula,
    Named("categories_order") = categories_order,
    Named("proportional") = proportional,
    Named("N_cats") = N_cats,
    Named("nobs_glmcat") = N,
    Named("distribution") = distribution,
    Named("freedom_degrees") = freedom_degrees
  );

  output_list.attr("class") = "glmcat";

  return output_list;
}


//' GLMcat model predictions
//'
//' @param model_object a GLMcat model
//' @param data a data frame in which to look for variables with which to predict_glmcat. Note that all predict_glmcator variables should be
//' present having the same names as the variables used to fit the model.
//' @param type he type of predict_glmcations. \code{"prob"} gives probabilities,
//' \code{"cum.prob"} gives cumulative probabilities and \code{"linear.predict"} gives
//' the linear predict_glmcator.
//' @rdname predict_glmcat
//' @export
//' @examples
//' data(DisturbedDreams)
//' mod1 <- GLMcat(formula = Level ~ Age,
//' data = DisturbedDreams, distribution = "logistic")
//' predict_glmcat(mod1, data = DisturbedDreams[1:5, ], type = "prob")
// [[Rcpp::export("predict_glmcat")]]
NumericVector predict_glmcat(List model_object,
                      DataFrame data,
                      String type
                      ){

  class distribution dist1;
  // Environment base_env("package:base");
  // Function my_rowSums = base_env["rowSums"];
  int N_cats = model_object["N_cats"];
  Eigen::MatrixXd coef = model_object["coefficients"];

  List newdataList = dist1.All_pre_data_NEWDATA(model_object["formula"],
                                                data,
                                                model_object["categories_order"],
                                                            model_object["proportional"],
                                                                        N_cats);

  Eigen::MatrixXd Design_Matrix = newdataList["Design_Matrix"];
  Eigen::MatrixXd predict_glmcated_eta;

  std::string ratio = model_object["ratio"];
  String distribution = model_object["distribution"];
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
      if(distribution == "logistic"){
        pi = ref.inverse_logistic(predict_glmcated_eta);
      }else if(distribution == "normal"){
        pi = ref.inverse_normal(predict_glmcated_eta);
      }else if(distribution == "cauchit"){
        pi = ref.inverse_cauchit(predict_glmcated_eta);
      }else if(distribution == "gompertz"){
        pi = ref.inverse_gompertz(predict_glmcated_eta);
      }else if(distribution == "gumbel"){
        pi = ref.inverse_gumbel(predict_glmcated_eta);
      }else if(distribution == "student"){
        pi = ref.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
      }
    }else if(ratio == "adjacent"){
      AdjacentR adj;
      if(distribution == "logistic"){
        pi = adj.inverse_logistic(predict_glmcated_eta);
      }else if(distribution == "normal"){
        pi = adj.inverse_normal(predict_glmcated_eta);
      }else if(distribution == "cauchit"){
        pi = adj.inverse_cauchit(predict_glmcated_eta);
      }else if(distribution == "gompertz"){
        pi = adj.inverse_gompertz(predict_glmcated_eta);
      }else if(distribution == "gumbel"){
        pi = adj.inverse_gumbel(predict_glmcated_eta);
      }else if(distribution == "student"){
        pi = adj.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
      }
    }else if(ratio == "sequential"){
      SequentialR seq;
      // Vector pi depends on selected distribution
      if(distribution == "logistic"){
        pi = seq.inverse_logistic(predict_glmcated_eta);
      }else if(distribution == "normal"){
        pi = seq.inverse_normal(predict_glmcated_eta);
      }else if(distribution == "cauchit"){
        pi = seq.inverse_cauchit(predict_glmcated_eta);
      }else if(distribution == "gompertz"){
        pi = seq.inverse_gompertz(predict_glmcated_eta);
      }else if(distribution == "gumbel"){
        pi = seq.inverse_gumbel(predict_glmcated_eta);
      }else if(distribution == "student"){
        pi = seq.inverse_student(predict_glmcated_eta, freedom_degrees);
      }else{
        Rcpp::stop("Unrecognized distribution; options are: logistic, normal, cauchit, gumbel, gompertz, and student(df)");
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

  if(type == "prob"){
    predict_glmcat = pi_total;
  }else if(type == "cum.prob"){
    predict_glmcat = cum_prob;
  }else if(type == "linear.predict"){
    predict_glmcat = Design_Matrix * coef;
  }else{
    Rcpp::stop("Unrecognized type parameter; options are: prob, cum.prob, linear.predict");
  }

  return predict_glmcat;
    // List::create(
    // Named("Design_Matrix") = Design_Matrix,
    // Named("Eta") = predict_glmcated_eta,
    // Named("cum_prob") = cum_prob,
    // Named("pi_total") = pi_total
    // Named("")
  // );

}


RCPP_MODULE(GLMcatmodule){
  Rcpp::function("GLMcat", &GLMcat,
                 List::create(_["formula"],
                              _["ratio"] = "reference",
                              _["distribution"] = "logistic",
                              _["categories_order"] = CharacterVector::create(NA_STRING),
                               _["proportional"] = CharacterVector::create(NA_STRING),
                               _["data"],
                                _["freedom_degrees"] = 1),
                                "GLMcat models");
  Rcpp::function("predict_glmcat", &predict_glmcat,
                 List::create(_["model_object"] = R_NaN,
                              _["data"],
                              _["type"] = "prob"
                 ),
                 "GLMcat model predictions");
}
