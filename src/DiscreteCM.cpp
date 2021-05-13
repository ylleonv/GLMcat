#include "cdf.h"
#include "reference.h"

using namespace std;
using namespace Rcpp ;
using namespace Eigen;

//' Family of models for Discrete Choice
//' @description Discrete choice model: Requires data in long form.
//' For each individual (or decision maker), there are multiple observations (rows),
//' one for each of the alternatives the individual could have chosen.
//' We call the group of observations for an individual a “case”.
//' Each case represents a single statistical observation although it comprises
//' multiple observations.
//' @title Discrete_CM
//' @rdname Discrete_CM
//' @name Discrete_CM
//' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model. A particularity for the formula is that for the case-specific variables, the user can define a specific effect for a category.
//' @param case_id a string with the name of the column that identifies each case.
//' @param alternatives a string with the name of the column that identifies the vector of alternatives the individual could have chosen.
//' @param reference a string indicating the reference category
//' @param alternative_specific a character vector with the name of the explanatory variables that are different for each case, these are the alternative specific variables. By default, the case specific variables are the explanatory variables that are not identify in here, but that are part of the formula.
//' @param data a dataframe (in a long format) object in R, with the dependent variable as factor.
//' @param cdf
//' \describe{
//' \item{\code{cdf}:}{a string indicating the F cdf, options are: logistic, normal, cauchy, student (any df), noncentralt, gompertz, gumbel and laplace.}
//' \item{\code{df}:}{an integer with the degrees of freedom of the 'cdf'}
//' \item{\code{mu}:}{an integer with the mu parameter of the 'cdf'}
//' }
//' @param intercept if "conditional" then the design will be equivalent to the conditional logit model
//' @param normalization the quantile to use for the normalization of the estimated coefficients where the logistic distribution is used as the base cumulative distribution function.
//' @param control
//' \describe{
//' \item{\code{maxit}:}{the maximum number of iterations for the Fisher scoring algorithm.}
//' \item{\code{epsilon}:}{a double with to fix the epsilon value}
//' \item{\code{beta_init}:}{an appropiate sized vector for the initial iteration of the algorithm}
//' }
//' @examples
//' library(GLMcat)
//' data(TravelChoice)
//' Discrete_CM(formula = choice ~ hinc + gc + invt,
//' case_id = "indv",alternatives = "mode", reference = "air",
//' data = TravelChoice,  alternative_specific = c("gc", "invt"),
//' cdf = "logistic")
//' @note For these models it is not allowed to exclude the intercept.
//' @export
// [[Rcpp::export("Discrete_CM")]]
List Discrete_CM(Formula formula,
                 String case_id,
                 String alternatives,
                 CharacterVector reference,
                 CharacterVector alternative_specific,
                 DataFrame data,
                 List cdf,
                 String intercept,
                 double normalization,
                 List control
){
  // std::string cdf2 = cdf[0];
  //
  // Rcout << cdf2 << std::endl;

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


  class cdf dist1;

  List Full_M = dist1.select_data_nested(formula,
                                         case_id,
                                         alternatives,
                                         reference,
                                         alternative_specific,
                                         data,
                                         intercept
                                           //   ,
                                           // ratio
  );

  Eigen::MatrixXd Y_init = Full_M["Response_M"] ;

  int ai = Y_init.row(0).sum();
  // Y_init = Y_init + MatrixXd::Ones(Y_init.rows(), Y_init.cols());

  if(!(ai == 0 | ai == 1)){
    Y_init = Y_init + MatrixXd::Ones(Y_init.rows(), Y_init.cols());
  }

  Rcout << "Y_init" << std::endl;
  // Rcout << ai << std::endl;
  Rcout << Y_init << std::endl;
  //
  // Rcout << "Y_init" << std::endl;

  Eigen::MatrixXd X_EXT = Full_M["Design_Matrix"];

  CharacterVector Names_design = Full_M["Names_design"];

  CharacterVector categories_order = Full_M["categories_order"];

  int Q = Y_init.cols();
  int K = Q + 1;
  int N = K * Y_init.rows();

  // Eigen::MatrixXd BETA;
  // BETA = Eigen::MatrixXd::Zero(X_EXT.cols(),1); // Beta initialization with zeros

  Eigen::MatrixXd BETA2;
  BETA2 = Eigen::MatrixXd::Zero(X_EXT.cols(),1);
  Eigen::VectorXd BETA3 = BETA2;

  Eigen::MatrixXd BETA = BETA2;

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

  double epsilon = 0.0001 ;
  double qp , s0 = 1;


  NumericVector beta_init;
  int iterations_us = 25;

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
  // Rcout << N/K << std::endl;

  if(beta_init.length() >= 2 ){
    // BETA = beta_init;
    BETA = as<Eigen::Map<Eigen::VectorXd> >(beta_init);
  }

  while ((Stop_criteria >( epsilon / N) ) & (iteration < ( iterations_us )) ){

    Eigen::MatrixXd Score_i = Eigen::MatrixXd::Zero(BETA.rows(),1);
    Eigen::MatrixXd F_i = Eigen::MatrixXd::Zero(BETA.rows(), BETA.rows());
    LogLik = 0.;

    // Rcout << cdf_1 << std::endl;
    // Rcout << Score_i << std::endl;
    // Rcout << F_i << std::endl;


    for (int i=0; i < N/K; i++){



      X_M_i = X_EXT.block(i*Q , 0 , Q , X_EXT.cols());
      // Rcout << "X_M_i" << std::endl;
      // Rcout << X_M_i << std::endl;
      Y_M_i = Y_init.row(i);

      // Rcout << "Y_M_i" << std::endl;
      // Rcout << Y_M_i << std::endl;

      eta = X_M_i * BETA;

      ReferenceF ref;
      if(cdf_1 == "logistic"){

        // Rcout << eta << std::endl;
        // print("as");
        pi = ref.inverse_logistic(eta);
        D = ref.inverse_derivative_logistic(eta);

        // Rcout << pi << std::endl;


      }else if(cdf_1 == "normal"){

        // Rcout << eta << std::endl;

        pi = ref.inverse_normal(eta);
        D = ref.inverse_derivative_normal(eta);

        // Rcout << pi << std::endl;

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
        Rcpp::stop("Unrecognized cdf; options are: logistic, normal, cauchy, gumbel, gompertz, laplace, student(df), and noncentral(df,mu)");
      }
      // }else{
      //   AdjacentR adj;
      //   if(cdf_1 == "logistic"){
      //     pi = adj.inverse_logistic(eta);
      //     D = adj.inverse_derivative_logistic(eta);
      //   }else if(cdf_1 == "normal"){
      //     pi = adj.inverse_normal(eta);
      //     D = adj.inverse_derivative_normal(eta);
      //   }else if(cdf_1 == "cauchy"){
      //     pi = adj.inverse_cauchy(eta);
      //     D = adj.inverse_derivative_cauchy(eta);
      //   }else if(cdf_1 == "student"){
      //     pi = adj.inverse_student(eta, freedom_degrees);
      //     D = adj.inverse_derivative_student(eta, freedom_degrees);
      //   }
      // }

      Cov_i = Eigen::MatrixXd(pi.asDiagonal()) - (pi*pi.transpose());
      W_in = D * Cov_i.inverse();
      Score_i_2 = (X_M_i.transpose() * W_in) * (Y_M_i - pi);
      Score_i = Score_i + Score_i_2;
      F_i_2 = X_M_i.transpose() * (W_in) * (D.transpose() * X_M_i);
      F_i = F_i + F_i_2;
      LogLik = LogLik + (Y_M_i.transpose().eval()*Eigen::VectorXd(pi.array().log())) + ( (1 - Y_M_i.sum()) * std::log(1 - pi.sum()) );

    }
    // To stop when LogLik is smaller than the previous
    if(iteration>5){
      if (LogLikIter[iteration] > LogLik)
        break;
    }

    LogLikIter.conservativeResize( LogLikIter.rows() +1 , 1);
    LogLikIter(LogLikIter.rows() - 1) = LogLik;
    Stop_criteria = (abs(LogLikIter(iteration+1) - LogLikIter(iteration))) / (epsilon + (abs(LogLikIter(iteration+1)))) ;
    VectorXd beta_old = BETA;

    // MatrixXd inverse;
    FullPivLU<MatrixXd> lu(F_i);
    bool invertible = lu.isInvertible();

    if(!invertible) {
      Rcpp::stop("Fisher matrix is not invertible");
    }

    // Rcout << "BETA" << std::endl;
    // Rcout << BETA << std::endl;

    MatrixXd F_inv = F_i.inverse();
    F_inv = F_inv * Score_i;
    // Rcout << F_inv << std::endl;

    BETA = BETA + F_inv;

    iteration = iteration + 1;
    // Rcout << "BETA" << std::endl;
    // Rcout << BETA << std::endl;

    // Rcout << "iteration" << std::endl;
    // Rcout << iteration << std::endl;
    // Rcout << "LogLik" << std::endl;
    // Rcout << LogLik << std::endl;


    F_i_final = F_i;
    // Rcout << "BETA" << std::endl;
    // Rcout << BETA << std::endl;
    // Rcout << "LogLik" << std::endl;
    // Rcout << LogLik << std::endl;

  }

  cov_beta = F_i_final.inverse();
  Std_Error = cov_beta.diagonal();
  Std_Error = Std_Error.array().sqrt() ;

  // Eigen::MatrixXd X_M_i_1 = X_EXT.block(0*Q , 0 , Q , X_EXT.cols());
  // Eigen::VectorXd Y_M_i_1 = Y_init.row(0);

  NumericMatrix BETA_2 = wrap(BETA);
  rownames(BETA_2) = Names_design;



  // Rcout << normalization << std::endl;
  // Rcout << cdf_1 << std::endl;

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

    NumericMatrix BETA_3 = BETA_2 * (s0);
    // output_list_dis["normalized_coefficients"] = BETA_3;
  }
  int df = BETA_2.length();

  std::string Convergence = "False";

  if(iteration < iterations_us){
    Convergence = "True";
  }


  List output_list_dis = List::create(
    Named("Function") = "DiscreteCM",
    Named("formula") = formula,
    Named("convergence") = Convergence,
    Named("ratio") = "reference",
    Named("Nb. iterations") = iteration-1 ,
    Named("coefficients") = BETA_2,
    Named("LogLikelihood") = LogLikIter(LogLikIter.rows() - 1),
    Named("LogLikIter") =  LogLikIter,
    Rcpp::Named("df of the model") = df,
    Named("X_M_i") =  X_M_i,
    Named("stderr") =  Std_Error,
    Named("N_cats") = K,
    Named("normalization_s0") =  s0,
    Named("cdf") = cdf,
    Named("nobs_glmcat") = N/K,
    Named("control") = control,
    Named("arguments") = List::create(Named("formula")= formula,Named("case_id")= case_id, Named("alternatives") = alternatives,
          Named("reference") = reference, Named("alternative_specific") = alternative_specific, Named("intercept") = intercept,
                 Named("categories_order") = categories_order)
  );


  output_list_dis.attr("class") = "glmcat";
  return output_list_dis;
}

RCPP_MODULE(Discrete_CMmodule){
  Rcpp::function("Discrete_CM", &Discrete_CM,
                 List::create(_["formula"] = R_NaN,
                              _["case_id"] = "a",
                              _["alternatives"] = "a",
                              _["reference"] = R_NaN,
                              _["alternative_specific"] = CharacterVector::create( NA_STRING),
                              _["data"] = NumericVector::create( 1, NA_REAL, R_NaN, R_PosInf, R_NegInf),
                              _["cdf"] = R_NaN,
                              _["intercept"] = "standard",
                              _["normalization"] = 1.0,
                              _["control"] = R_NaN
                 ),
                 "Discrete Choice Model");

}
