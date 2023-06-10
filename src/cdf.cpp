
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]

#include "cdf.h"
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/non_central_t.hpp>


using namespace boost::math;
using namespace std;
using namespace Rcpp ;

cdf::cdf(void) {
  // Rcout << "cdf is being created" << endl;
}



// std::string cdf::concatenate(std::string x, std::string level)
// {
//   return (x + " " +level);
// }

// ORDINAL DATA SELECTION

Environment base_base("package:base");
Function my_transpose1 = base_base["t"];

std::string cdf::concatenate(std::string x, std::string level)
{
  return (x + " " +level);
}

// Environment stats_env1("package:utils");
// Function head1 = stats_env1["head"];

DataFrame my_transpose(DataFrame dat_in){
  // Environment base_base("package:base");
  // Function my_transpose1 = base_base["t"];
  DataFrame data_out = my_transpose1(dat_in);
  return data_out;
}

NumericMatrix to_dummy1(NumericVector A, CharacterVector levels)
{
  NumericVector cha_to_fact = A;
  CharacterVector levs1 = levels;
  int var_lev = levs1.length();

  // Rcout << var_lev << std::endl;

  NumericMatrix B_Ma(cha_to_fact.length(), var_lev);
  for (int i_1 = 0; i_1 < cha_to_fact.length(); ++i_1){
    int col_ind = cha_to_fact[i_1] - 1;
    B_Ma(i_1, col_ind) = 1;
  }

  // print(B_Ma);

  B_Ma = B_Ma( _ , Range(0,var_lev-2) );

  // if (var_lev != 2){
  //   B_Ma = B_Ma( _ , Range(0,var_lev-2) );
  // } else {
  //   B_Ma = B_Ma( _ , Range(1,var_lev-1) );
  // }
  return B_Ma;
}

// [[Rcpp::export("Model_Matrix_or")]]
List Model_Matrix_or(DataFrame data, Formula formula) {
  Environment stats_env("package:stats");
  Environment base_env("package:base");
  Function model_frame = stats_env["model.frame"];
  Function model_matrix = stats_env["model.matrix"];

  Function names_f = base_env["names"];
  DataFrame df_new1 = model_frame(Rcpp::_["formula"] = formula, _["data"] = data);
  SEXP Response = df_new1[0];
  CharacterVector covariates_names = names_f(df_new1);

  // Function n_levels = base_env["nlevels"];
  Function is_factor = base_env["is.factor"];
  Function is_integer = base_env["is.integer"];
  // Function as_numeric = base_env["as.numeric"];
  Function sapply_f = base_env["sapply"];
  // Function lapply_f = base_env["lapply"];


  // Variables which are integers are converted to numeric
  NumericVector int_var = sapply_f(data[covariates_names], is_integer);
  int_var = abs(int_var - 1);
  SEXP int_var1 = int_var;
  LogicalVector int_var2 = int_var1;
  CharacterVector var_int = covariates_names[int_var2];
  // Rcout << int_var << std::endl;
  // Rcout << var_int << std::endl;


  // CharacterVector covariates_names = M_matrix["covariates_names"];
  LogicalVector factor_var1 = sapply_f(data[var_int], is_factor);
  // factor_var = factor_var[var_int];
  CharacterVector factor_var = var_int[factor_var1];
  // Rcout << "factor_var" << std::endl;
  // Rcout << df_new << std::endl;

  NumericMatrix df_new = model_matrix(df_new1, _["data"] = df_new1);

  // Rcout << df_new << std::endl;

  return List::create(
    Named("df_new") = df_new,
    Named("Response") = Response,
    Named("covariates_names") = covariates_names,
    Named("factor_var") = factor_var
  );
}

List Model_Matrix_or_pre(DataFrame data, Formula formula) {
  // modification of the formula and just model.matrix with drop
  Environment stats_env("package:stats");
  Environment base_stats("package:base");
  Function model_frame = stats_env["model.frame"];
  Function model_matrix = stats_env["model.matrix"];
  Function my_strsplit = base_stats["strsplit"];
  Function my_paste = base_stats["paste"];
  Function my_as_formula = stats_env["as.formula"];
  Function my_format = base_stats["format"];
  // Here I take off the response to have only the explanatory variables as entry to the formula ~ predictors
  // in order to predict for new data
  CharacterVector st1 = my_format(formula);
  String str_for = my_paste(st1,  _["collapse"] = "");
  List list1 = (my_strsplit(str_for, "~"));
  CharacterVector vars = list1[0];
  String vars_string = vars[1];
  Formula new_for = my_as_formula(my_paste("~", vars_string,  _["collapse"] = ""));

  DataFrame df_new1 = model_frame(Rcpp::_["formula"] = new_for, _["data"] = data);
  SEXP Response = df_new1[0];
  NumericMatrix df_new = model_matrix(df_new1, _["data"] = data);
  return List::create(
    Named("df_new") = df_new,
    Named("Response") = Response
  );
}


List Cat_ref_order(CharacterVector categories_order, SEXP response_categories){
  // Environment base_env("package:base");
  // Function my_unique = base_env["unique"];

  CharacterVector response_categories1 = response_categories;
  // CharacterVector levels = my_unique(response_categories1);

  IntegerVector num_categories_order = seq_len(categories_order.length());

  DataFrame response_neworder = DataFrame::create(  _["num_categories_order"] = num_categories_order );

  response_neworder = my_transpose(response_neworder);
  response_neworder.names() = categories_order;

  String a1;
  CharacterVector a2;
  for(int i = 0; i <response_categories1.length(); i++){
    a1 = response_categories1[i];
    a2 = response_neworder[a1];
    response_categories1[i] = a2[0];
  }

  // print(response_categories1);

  return List::create(
    Named("response_neworder") = response_neworder,
    Named("response_categories2") = response_categories1,
    Named("levels") = categories_order
  );
}



List cdf::All_pre_data_or(Formula formula,
                          DataFrame input_data,
                          CharacterVector categories_order,
                          CharacterVector parallel_effect,
                          std::string threshold,
                          std::string ratio){

  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_cbind = base_env["cbind"];
  Function my_order = base_env["order"];
  Function my_levels = base_env["levels"];

  List M_matrix = Model_Matrix_or(input_data, formula);

  if(categories_order.length() == 1){
    LogicalVector is_na_ref = is_na(categories_order);
    // Rcout << is_na_ref << std::endl;
    // print(my_levels(M_matrix["Response"]));

    if(!is_na_ref[0]){ // If it is a character diferent of NA
      CharacterVector categories_order_n = my_levels(M_matrix["Response"]);
      // Rcout << categories_order_n << std::endl;
      int ind_ref = 0;
      while(categories_order[0] != categories_order_n(ind_ref)){ind_ref = ind_ref+1;}
      categories_order_n.erase(ind_ref);
      categories_order_n.push_back(categories_order[0]);
      // Rcout << categories_order_n << std::endl;
      categories_order = categories_order_n;
    }else{ // For other case default alphabetical
      categories_order = my_levels(M_matrix["Response"]);
    }
  }

  // Rcout << categories_order_n << std::endl;

  List Cat_ref_or_L = Cat_ref_order(categories_order, M_matrix["Response"]);
  NumericMatrix Design = M_matrix["df_new"];

  CharacterVector a1 = Cat_ref_or_L["response_categories2"];
  NumericVector Num_res = my_asnumeric(a1);
  Design = my_cbind(Num_res, Design);

  // Rcout << "is_na_ref1" << std::endl;
  // Now order dataset with respect to the repsonse variables in the given order
  DataFrame df_tans = my_transpose(Design);
  DataFrame df_tans_2 = Design ;
  // NumericVector order_var_sel = my_order(df_tans_2[0]);
  // order_var_sel = order_var_sel - 1 ;
  // df_tans = df_tans[order_var_sel];
  df_tans_2 = my_transpose(df_tans);
  CharacterVector Levels = Cat_ref_or_L["levels"];
  int N_cats = Levels.length();

  // Rcout << "is_na_ref2" << std::endl;


  CharacterVector colnames_final_m = df_tans_2.names();
  CharacterVector factor_var = M_matrix["factor_var"];
  CharacterVector covariates_names = M_matrix["covariates_names"];

  LogicalVector any_alternative_specific = is_na(parallel_effect); // TRUE IF NA

  // Rcout << any_alternative_specific << std::endl;
  // Rcout << parallel_effect << std::endl;

  // Assume parallel design
  if((parallel_effect[0] == "TRUE" || parallel_effect[0] == "NA") && (parallel_effect.size() == 1)){
    // Rcout << parallel_effect.size() << std::endl;
    // if(parallel_effect[0] == "TRUE" && (parallel_effect.size() ==1)){
    parallel_effect = colnames_final_m;
    parallel_effect.erase(0);
    parallel_effect.erase(0);
    // Rcout << parallel_effect << std::endl;
  }else if(parallel_effect[0] == "FALSE" && (parallel_effect.size() ==1)){ // Explicit complete design
    parallel_effect = CharacterVector::create(NA_STRING);
    // Rcout << parallel_effect << std::endl;
  }

  // Rcout << factor_var << std::endl;

  any_alternative_specific = !is_na(parallel_effect); // TRUE IF THERE ARE

  CharacterVector parallel_effect1 = parallel_effect;
  int con = 0;

  // Identify which parallel effects are factors
  NumericVector x1(parallel_effect.length());
  if (any_alternative_specific[0]) {
    for(int indi = 0 ; indi < parallel_effect.length(); indi++){
      String var_1 = parallel_effect[indi];
      for(int id_fac = 0; id_fac < factor_var.length(); id_fac++){
        if(factor_var[id_fac] == var_1){
          x1[indi]=1;
          parallel_effect1.erase(indi-con);
          con = con + 1;
        }
      }
    }
  }

  SEXP x2 = x1;
  LogicalVector x3 = x2;
  CharacterVector par_factor = parallel_effect[x3];
  // print(x1);
  // Rcout << par_factor << std::endl;

  Function f_levels = base_env["levels"];
  Function f_c = base_env["c"];
  Function my_paste = base_base["paste"];

  CharacterVector names_factors; // Extended names of parallel covariates that are factors

  for(int indi = 0 ; indi < par_factor.length(); indi++){
    String var_1 = par_factor[indi];
    CharacterVector levels_fac = f_levels(input_data[var_1]);
    levels_fac.erase(0); // Reference is the first level
    CharacterVector names_factor = my_paste(var_1, levels_fac, _["sep"] = "");
    names_factors = f_c(names_factors, names_factor);
  }
  // Rcout << names_factors << std::endl;

  parallel_effect = f_c(parallel_effect1,names_factors);
  // Rcout << parallel_effect << std::endl;

  if (any_alternative_specific[0]) {
  for(int indi = 0 ; indi < parallel_effect.length(); indi++){
    String var_1 = parallel_effect[indi];
    string victor = var_1;
    const char* ccx = victor.c_str();
    if(df_tans_2.containsElementNamed(ccx)){
    }else{
      parallel_effect.erase(indi);
    }
  }}

  // Rcout << parallel_effect << std::endl;

  NumericVector x(colnames_final_m.length());
  if (any_alternative_specific[0]) {
    // x(colnames_final_m.length());
    for(int indi = 0 ; indi < parallel_effect.length(); indi++){
      String var_1 = parallel_effect[indi];
      // if(df_tans_2.containsElementNamed(ccx)){
      int indi_var = df_tans_2.findName(var_1);
      // Rcout << "is_na_ref4" << std::endl;
      x[indi_var] = indi_var;
      // }else{
      // parallel_effect.erase(indi_var);
      // Rcout << "is_na_ref3" << std::endl;
      // }
    }
    // print(colnames_final_m);
    colnames_final_m = colnames_final_m[x==0]; // Case where there no parallel effects
    // print(colnames_final_m);
  }

  // Rcout << "is_na_ref4" << std::endl;

  // Now extend
  // Y EXTEND
  NumericVector Response = df_tans_2[0];
  NumericMatrix Response_EXT = to_dummy1(Response, Levels);


  DataFrame DF_complete_effect = df_tans_2[colnames_final_m];
  NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  Eigen::MatrixXd X_EXT_COMPLETE;
  // Rcout << Pre_Design << std::endl;
  // print(head1(Pre_Design));

  if (ratio == "cumulative"){
    if (threshold == "equidistant"){
      // if (!any_alternative_specific[0]) { // ninguna es parallel
      // Rcout << "propor_cum_equ" << std::endl;

      NumericMatrix tJac = my_cbind(1, seq_len(categories_order.length() -1 )-1 );
      Eigen::Map<Eigen::MatrixXd> tJac2 = as<Eigen::Map<Eigen::MatrixXd> >(tJac);
      // Rcout << tJac2 << std::endl;
      Eigen::MatrixXd X_EXT_COMPLETE_int = Eigen::kroneckerProduct(Pre_Design.block(0,1,Pre_Design.rows(),1), tJac2).eval();
      Eigen::MatrixXd X_EXT_COMPLETE_res = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 2), Iden_Q1).eval();

      Eigen::MatrixXd X_EXT_COMPLETE_1 = Eigen::MatrixXd::Zero(X_EXT_COMPLETE_int.rows(),
                                                               X_EXT_COMPLETE_int.cols()+X_EXT_COMPLETE_res.cols());

      X_EXT_COMPLETE_1.block(0,0,X_EXT_COMPLETE_int.rows(),X_EXT_COMPLETE_int.cols()) = X_EXT_COMPLETE_int;
      X_EXT_COMPLETE_1.block(0,X_EXT_COMPLETE_int.cols(),X_EXT_COMPLETE_res.rows(),X_EXT_COMPLETE_res.cols()) = X_EXT_COMPLETE_res;

      // print(head1(X_EXT_COMPLETE_1));
      X_EXT_COMPLETE = X_EXT_COMPLETE_1;

      colnames_final_m.erase(0);
      // Rcout << X_EXT_COMPLETE << std::endl;
      // }
    }else if (threshold == "symmetric"){ // caso symmetric

      Eigen::MatrixXd Sim_int = Eigen::MatrixXd::Ones(N_cats-1,floor(N_cats/2)+1);
      Eigen::MatrixXd Iden_sim = -Eigen::MatrixXd::Identity(floor(N_cats/2),floor(N_cats/2)).rowwise().reverse();
      // Rcout << Iden_sim  << std::endl;

      Eigen::MatrixXd Iden_sim1 = Eigen::MatrixXd::Identity(floor(N_cats/2),floor(N_cats/2));
      // Rcout << Iden_sim1  << std::endl;

      Sim_int.block(0,1,Iden_sim.rows(),Iden_sim.cols()) = Iden_sim;
      Sim_int.block(N_cats-floor(N_cats/2)-1,1,Iden_sim1.rows(),Iden_sim1.cols()) = Iden_sim1;

      if((N_cats % 2) == 0) { //for even numbers the row in the middle to zero execpt the intercept
        Sim_int.block((N_cats/2) - 1 ,1,1,Sim_int.cols()-1) = Eigen::VectorXd::Zero(Sim_int.cols()-1);
      }

      // Rcout << Sim_int  << std::endl;

      Eigen::MatrixXd X_EXT_COMPLETE_int = Eigen::kroneckerProduct(Pre_Design.block(0,1,Pre_Design.rows(),1), Sim_int).eval();
      Eigen::MatrixXd X_EXT_COMPLETE_res = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 2), Iden_Q1).eval();

      Eigen::MatrixXd X_EXT_COMPLETE_1 = Eigen::MatrixXd::Zero(X_EXT_COMPLETE_int.rows(),
                                                               X_EXT_COMPLETE_int.cols()+X_EXT_COMPLETE_res.cols());

      X_EXT_COMPLETE_1.block(0,0,X_EXT_COMPLETE_int.rows(),X_EXT_COMPLETE_int.cols()) = X_EXT_COMPLETE_int;
      X_EXT_COMPLETE_1.block(0,X_EXT_COMPLETE_int.cols(),X_EXT_COMPLETE_res.rows(),X_EXT_COMPLETE_res.cols()) = X_EXT_COMPLETE_res;

      X_EXT_COMPLETE = X_EXT_COMPLETE_1;

      colnames_final_m.erase(0);

    }else{ //standard
      X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();
      colnames_final_m.erase(0);
    }
  }else{ // caso not cum
    X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();
    colnames_final_m.erase(0);
  }

  Eigen::MatrixXd Design_Matrix;


  // print(parallel_effect);

  if (any_alternative_specific[0]) { // para las proporcionales
    DataFrame DF_parallel_effect = df_tans_2[parallel_effect];
    NumericMatrix Pre_DF_parallel1 = internal::convert_using_rfunction(DF_parallel_effect, "as.matrix");
    Eigen::Map<Eigen::MatrixXd> Pre_DF_parallel2 = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_parallel1);
    Eigen::MatrixXd Pre_DF_parallel = Pre_DF_parallel2;
    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
    Eigen::MatrixXd X_EXT_parallel = Eigen::kroneckerProduct(Pre_DF_parallel, Ones).eval();

    Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_parallel.cols());
    Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
    Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_parallel.cols()) = X_EXT_parallel;
  }else{
    Design_Matrix = X_EXT_COMPLETE;
  }
  // Rcout << Design_Matrix << std::endl;
  // print(colnames_final_m);

  return List::create(
    Named("Design_Matrix") = Design_Matrix,
    Named("Response_EXT") = Response_EXT,
    Named("Levels") = Levels,
    Named("Complete_effects") = colnames_final_m,
    Named("parallel_effect") = parallel_effect,
    Named("N_cats") = N_cats,
    Named("categories_order") = categories_order
  );
}

List cdf::All_pre_data_NEWDATA(Formula formula,
                               DataFrame NEWDATA,
                               CharacterVector categories_order,
                               CharacterVector parallel_effect,
                               int N_cats){

  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_cbind = base_env["cbind"];
  Function my_order = base_env["order"];
  List M_matrix = Model_Matrix_or_pre(NEWDATA, formula);
  // List Cat_ref_or_L = Cat_ref_order(categories_order, M_matrix["Response"]);
  NumericMatrix Design = M_matrix["df_new"];

  // CharacterVector a1 = Cat_ref_or_L["response_categories2"];
  // NumericVector Num_res = my_asnumeric(a1);
  //
  // Design = my_cbind(Num_res, Design);

  // Now order dataset with respect to the repsonse variables in the given order
  // DataFrame df_tans = my_transpose(Design);
  DataFrame df_tans_2 = Design ;


  // NumericVector order_var_sel = my_order(df_tans_2[0]);
  // order_var_sel = order_var_sel - 1 ;
  // df_tans = df_tans[order_var_sel];

  // df_tans_2 = my_transpose(df_tans);

  // CharacterVector Levels = Cat_ref_or_L["levels"];
  // int N_cats = Levels.length();

  LogicalVector any_alternative_specific = !is_na(parallel_effect); // TRUE IF THERE ARE
  CharacterVector colnames_final_m = df_tans_2.names();


  NumericVector x(colnames_final_m.length());
  if (any_alternative_specific[0]) {
    // x(colnames_final_m.length());
    for(int indi = 0 ; indi < parallel_effect.length(); indi++){
      String var_1 = parallel_effect[indi];
      int indi_var = df_tans_2.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no parallel effects
  }


  // X EXTEND

  // X COMPLETE

  DataFrame DF_complete_effect = df_tans_2[colnames_final_m];

  NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);

  Eigen::MatrixXd X_EXT_COMPLETE;

  // if (threshold == "equidistant"){ // eRASE THE LAST TWO COLUMNS ONE CORRESPONDING TO DR AND OTHER TO THE INTERCEPT
  //   X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 2), Iden_Q1).eval();
  //   colnames_final_m.erase(0);
  //   colnames_final_m.erase(0);
  // } else {
  X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design, Iden_Q1).eval();


  // colnames_final_m.erase(0);
  // }

  Eigen::MatrixXd Design_Matrix;

  if (any_alternative_specific[0]) {

    // PONER ESA PARTE ACA

    DataFrame DF_parallel_effect = df_tans_2[parallel_effect];
    NumericMatrix Pre_DF_parallel1 = internal::convert_using_rfunction(DF_parallel_effect, "as.matrix");
    Eigen::Map<Eigen::MatrixXd> Pre_DF_parallel2 = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_parallel1);

    Eigen::MatrixXd Pre_DF_parallel = Pre_DF_parallel2;

    Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
    Eigen::MatrixXd X_EXT_parallel = Eigen::kroneckerProduct(Pre_DF_parallel, Ones).eval();

    // TENGO QUE PONER EL IF ACA
    //
    // if (threshold == "equidistant"){
    //   NumericMatrix tJac = my_cbind(1, seq_len(categories_order.length() -1 )-1 );
    //   Eigen::Map<Eigen::MatrixXd> tJac2 = as<Eigen::Map<Eigen::MatrixXd> >(tJac);
    //   Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(Response_EXT.rows());
    //   Eigen::MatrixXd Threshold_M = Eigen::kroneckerProduct(Ones1, tJac2).eval();
    //   X_EXT_parallel.conservativeResize(X_EXT_parallel.rows(),X_EXT_parallel.cols()+2);
    //   X_EXT_parallel.block(0,X_EXT_parallel.cols()-2, X_EXT_parallel.rows(),2) = Threshold_M;
    // }


    Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_parallel.cols());
    Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
    Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_parallel.cols()) = X_EXT_parallel;

  }else{Design_Matrix = X_EXT_COMPLETE;}

  // Now extend


  // // X COMPLETE
  // DataFrame DF_complete_effect = df_tans_2[colnames_final_m];
  // NumericMatrix Pre_Design1 = internal::convert_using_rfunction(DF_complete_effect, "as.matrix");
  // Eigen::Map<Eigen::MatrixXd> Pre_Design = as<Eigen::Map<Eigen::MatrixXd> >(Pre_Design1);
  // // print(colnames_final_m);
  // Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  //
  // Eigen::MatrixXd X_EXT_COMPLETE = Eigen::kroneckerProduct(Pre_Design.rightCols(DF_complete_effect.cols() - 1), Iden_Q1).eval();
  //
  // // X PROPOTIONAL
  //
  // Eigen::MatrixXd Design_Matrix;
  //
  // if (any_alternative_specific[0]) {
  //
  //   DataFrame DF_parallel_effect = df_tans_2[parallel_effect];
  //   NumericMatrix Pre_DF_parallel1 = internal::convert_using_rfunction(DF_parallel_effect, "as.matrix");
  //   Eigen::Map<Eigen::MatrixXd> Pre_DF_parallel = as<Eigen::Map<Eigen::MatrixXd> >(Pre_DF_parallel1);
  //   Eigen::VectorXd Ones = Eigen::VectorXd::Ones(N_cats-1);
  //   Eigen::MatrixXd X_EXT_parallel = Eigen::kroneckerProduct(Pre_DF_parallel, Ones).eval();
  //
  //   Design_Matrix.conservativeResize(X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()+X_EXT_parallel.cols());
  //   Design_Matrix.block(0,0,X_EXT_COMPLETE.rows(),X_EXT_COMPLETE.cols()) = X_EXT_COMPLETE;
  //   Design_Matrix.block(0,X_EXT_COMPLETE.cols(),X_EXT_COMPLETE.rows(),X_EXT_parallel.cols()) = X_EXT_parallel;
  //
  // }else{Design_Matrix = X_EXT_COMPLETE;}

  return List::create(
    Named("Design") = Design,
    Named("Design_Matrix") = Design_Matrix
  );
}


CharacterVector Var_Not_In(DataFrame final_matrix, CharacterVector alternative_specific){

  LogicalVector any_alternative_specific = !is_na(alternative_specific);

  CharacterVector colnames_final_m = final_matrix.names();

  if (any_alternative_specific[0]) {
    NumericVector x(colnames_final_m.length());
    for(int indi = 0 ; indi < alternative_specific.length(); indi++){
      String var_1 = alternative_specific[indi];
      int indi_var = final_matrix.findName(var_1);
      x[indi_var] = indi_var;
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no alternative specific variables
  }
  colnames_final_m.erase(0, 4); // Eliminate the first 4 information variables

  return colnames_final_m;
}

List formula_entry(Formula formula1){
  Environment base_base("package:base");
  Function my_strsplit = base_base["strsplit"];
  Function my_format = base_base["format"];
  Function my_paste = base_base["paste"];
  Function my_sub = base_base["sub"];
  Function my_rev = base_base["rev"];
  Function my_trimws = base_base["trimws"];

  Environment base_stats("package:stats");
  Function my_as_formula = base_stats["as.formula"];

  CharacterVector st1 = my_format(formula1);
  // print(st1);
  String str_for = my_paste(st1,  _["collapse"] = "");
  List list1 = (my_strsplit(str_for, "~"));
  CharacterVector vars = list1[0];
  String vars_string = vars[1];
  String res1 = vars[0];
  String res = my_trimws(res1);
  List list_vars = my_strsplit(vars_string, "[+]");
  CharacterVector char_vars = list_vars(0);
  CharacterVector char_vars1 = (my_trimws(char_vars));
  // LogicalVector intercept_logic = (my_trimws(char_vars) == "-1");
  // print(char_vars);

  LogicalVector intercept_logic(char_vars.size());
  for( int i=0; i<char_vars.size(); i++){
    intercept_logic[i] = ((char_vars1[i] != "-1") && (char_vars1[i] != "1" ));
  }
  bool yes_intercept = is_true(all(intercept_logic)); // ARE ALL TRUE. NEITHER EQUALS TO -1 or 1
  if(yes_intercept){
    char_vars.push_front("1");
  }

  String vars_for = my_paste(my_sub("\\[[^()]*\\]", "", char_vars),  _["collapse"] = " + ");
  CharacterVector form = my_paste(res, vars_for, _["sep"] = "~");
  Formula formula2 = my_as_formula(form);

  String firs_col;
  List list_no_spaces, list_var_rev, list_cat_rev;
  StringVector Vars(char_vars.length()), Alternatives(char_vars.length());

  for (int i = 0; i < char_vars.length() ; i++) {
    firs_col = char_vars[i];
    String no_spaces = my_trimws(firs_col);
    list_no_spaces = my_strsplit(no_spaces, "");
    String rev_string = my_paste(my_rev(list_no_spaces[0]), _["collapse"] = "");

    String var = my_sub(".*\\[", "", rev_string);
    list_var_rev = my_strsplit(var, "");
    String var1 = my_paste(my_rev(list_var_rev[0]), _["collapse"] = "");

    String cat = my_sub(".*\\]", "", my_sub("\\[.*", "", rev_string));
    list_cat_rev = my_strsplit(cat, "");
    String cat1 = my_paste(my_rev(list_cat_rev[0]), _["collapse"] = "");

    if (var1 != cat1){ Alternatives[i] = cat1; }
    if (var1 == "1"){ var1 = "(Intercept)"; }
    Vars[i] = var1;
  }
  DataFrame Var_alt = DataFrame::create(Named("Alternatives") = Alternatives);
  Var_alt = my_transpose(Var_alt);
  Var_alt.names() = Vars;

  DataFrame Var_alt1 = Var_alt;

  DataFrame Var_alt2 = DataFrame::create(Named("Alternatives") = Alternatives, Named("Vars") = Vars);

  // Rcout << formula2 << std::endl;

  return List::create(
    Named("Var_alt") = Var_alt,
    Named("Response") = res,
    Named("formula_model") = formula2,
    Named("Var_alt2") = Var_alt2
  );
}


List Cat_ref1(CharacterVector categories_order,
              RObject response_categories){
  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];

  CharacterVector Levels1 = response_categories.attr("levels");
  // print(Levels1);
  CharacterVector response_categories1 = as<CharacterVector>(response_categories);
  // If categories_order only has one level then I add it to the others ordered alphabetically
  if(categories_order.length() == 1){
    String categories_order1 = categories_order[0];
    for(int var = 0 ; var < Levels1.length() ; var++){
      if (categories_order1 == Levels1[var]){Levels1.erase(var);}
    }
    Levels1.push_back(categories_order1);
    categories_order = Levels1;
  }


  IntegerVector num_categories_order = seq_len(categories_order.length());
  DataFrame response_neworder = DataFrame::create(  _["num_categories_order"] = num_categories_order );

  response_neworder = my_transpose(response_neworder);
  response_neworder.names() = categories_order;

  String a1;
  CharacterVector a2;
  for(int i = 0; i <response_categories1.length(); i++){
    a1 = response_categories1[i];
    a2 = response_neworder[a1];
    response_categories1[i] = a2[0];
  }

  // print(response_categories1);

  return List::create(
    Named("response_neworder") = response_neworder,
    Named("response_categories2") = my_asnumeric(response_categories1),
    Named("levels") = categories_order
  );
}

List Sort_DataFrame(DataFrame ModelMatrix,
                    DataFrame InputData,
                    CharacterVector names,
                    String Choice_vec,
                    CharacterVector Ref_cat,
                    String predict) {
  // Order DATAFRAME ACCORDING TO VARIABLES GIVEN IN VECTOR NAMES
  // CBIND OF DATA SETS AND THEN ORDER ACCORDING TO VARIABLES
  Environment base_env("package:base");
  Function my_order = base_env["order"];
  Function my_cbind = base_env["cbind"];
  Function my_asnumeric = base_env["as.numeric"];

  String alt = names[0];
  String id_case_0 = names[2];

  List Cat_ref_vec1 = Cat_ref1(Ref_cat, InputData[alt]);
  CharacterVector Cat_ref_vec = Cat_ref_vec1["response_categories2"];
  CharacterVector Levels = Cat_ref_vec1["levels"];

  DataFrame A2 = my_cbind( _["alternatives"] = InputData[alt],
                           _["Cat_ref_vec"] = Cat_ref_vec,
                           _["id_case"] = InputData[id_case_0],
                                                   _["choice"] = my_asnumeric(InputData[Choice_vec]),
                                                   ModelMatrix);

  CharacterVector names1 = CharacterVector::create("alternatives", "Cat_ref_vec", "id_case");

  DataFrame df_tans = my_transpose(A2);
  DataFrame df_tans_2 = A2 ;

  for (int i = 0; i < names1.length() ; i++) {
    String var = names1(i);
    NumericVector order_var_sel = my_order(df_tans_2[var]);
    order_var_sel = order_var_sel - 1 ;
    // # para el predict
    if(predict != "predict"){
      df_tans = df_tans[order_var_sel];
    }
    df_tans_2 = my_transpose(df_tans);
  }

  return List::create(
    _["df_tans_2"] = df_tans_2,
    _["Levels"] = Levels
  );
}


DataFrame my_AsNumericMatrix(DataFrame dat_in){
  Environment base_env("package:base");
  Function my_asnumeric = base_env["as.numeric"];
  Function my_ascharacter = base_env["as.character"];
  Function my_cbind = base_env["cbind"];
  DataFrame data_out = dat_in ;

  for (int i = 4; i < dat_in.length() ; i++) {
    NumericVector vec = my_asnumeric(my_ascharacter(dat_in[i]));
    data_out[i] = vec;
  }
  return data_out;
}

NumericMatrix Model_Matrix(DataFrame data, Formula formula) {
  Environment stats_env("package:stats");
  Function model_frame = stats_env["model.frame"];
  Function model_matrix = stats_env["model.matrix"];
  DataFrame df_new1 = model_frame(Rcpp::_["formula"] = formula, _["data"] = data);
  NumericMatrix df_new = model_matrix(df_new1, _["data"] = data);
  return  df_new;
}

// // [[Rcpp::export]]
List All_pre_data(Formula formula,
                  DataFrame input_data,
                  CharacterVector var_informatives,
                  String choice,
                  CharacterVector Ref_cat,
                  String predict){

  // print(formula_entry(formula)["formula_model"]);
  // print(head1(Model_Matrix(input_data, formula_entry(formula)["formula_model"])));

  List Out_SM = Sort_DataFrame(
    Model_Matrix(input_data, formula_entry(formula)["formula_model"]),
    input_data,
    var_informatives,
    choice,
    Ref_cat,
    predict);

  SEXP MA1 = Out_SM["df_tans_2"];
  CharacterVector Levels = Out_SM["Levels"];
  DataFrame MA11 = Rcpp::as<DataFrame>(MA1);
  DataFrame data_output = my_AsNumericMatrix(MA11);

  // print(head1(data_output));

  return List::create(
    _["data_output"] = data_output,
    _["Levels"] = Levels
  );
}

Eigen::MatrixXd Extend_alt_specific(DataFrame alt_specific, int N_cats, int N_ind,
                                    CharacterVector var_alt_specific
                                      // , String ad_or_ref
){
  Eigen::VectorXd Ones1 = Eigen::VectorXd::Ones(N_cats-1);
  Eigen::MatrixXd Iden_Q = Eigen::MatrixXd::Identity(N_cats,N_cats);

  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats+1,N_cats+1);
  Eigen::MatrixXd Iden_Q11 =  Iden_Q1.block(1, 0, N_cats, N_cats);
  // if(ad_or_ref != "reference"){Iden_Q = Iden_Q - Iden_Q11;
  //   Iden_Q.conservativeResize(Iden_Q.rows() - 1, Iden_Q.cols());
  // }else{
  Iden_Q.conservativeResize(Iden_Q.rows() - 1, Iden_Q.cols());
  Iden_Q.col(N_cats-1) = -Ones1;
  // }
  NumericMatrix alt_specific_num = internal::convert_using_rfunction(alt_specific[var_alt_specific], "as.matrix");
  Eigen::Map<Eigen::MatrixXd> M_alt_specific = as<Eigen::Map<Eigen::MatrixXd> >(alt_specific_num);
  Eigen::MatrixXd Matrix_trans((N_cats-1)*N_ind,var_alt_specific.length());
  for(int indi = 1 ; indi <= N_ind ; indi++)
  {
    Eigen::MatrixXd Block_ind =  M_alt_specific.block((indi-1) * N_cats, 0, N_cats, var_alt_specific.length());
    Eigen::MatrixXd Block_RES = Block_ind.transpose() * Iden_Q.transpose();
    Matrix_trans.block((indi-1) * (N_cats-1), 0, N_cats-1, var_alt_specific.length()) = Block_RES.transpose();
  }

  return Matrix_trans;
}

// Para cada individuo, como sus ingresos, seran los mismos para todas las categorias
Eigen::MatrixXd Extend_case_specific(DataFrame case_specific, int N_cats, int N_ind,
                                     CharacterVector var_alt_specific,  DataFrame Var_alt,
                                     String ref_cat, String intercept
                                       // , String ad_or_ref
){

  CharacterVector var_case_specific = Var_Not_In(case_specific, var_alt_specific);
  DataFrame effect_specific_for1 = Var_alt[var_case_specific];
  effect_specific_for1 = my_transpose(effect_specific_for1);
  CharacterVector effect_specific_for2 = effect_specific_for1[0];

  // 1 SI LA CATEGORIA DE REFERENCIA ES LA MISMA PARA LA QUE SE QUIERE EL EFECTO PARTICULAR EN LAS CASE-SPECIFIC
  NumericVector is_ref_alt_spe(effect_specific_for2.length());
  for (int i = 0; i < effect_specific_for2.length(); i++){
    if(ref_cat == effect_specific_for2(i)){
      is_ref_alt_spe(i) = 1;
    }else{is_ref_alt_spe(i) = 0;}
  }

  // Para hinc[bus] que tienen ref car
  Eigen::MatrixXd Iden_Q1 = Eigen::MatrixXd::Identity(N_cats-1,N_cats-1);
  // Para hinc[car] que tienen ref car
  Eigen::MatrixXd ONES_Q1 = Eigen::MatrixXd::Constant(N_cats-1,N_cats-1, -1);

  NumericMatrix case_specific_num = internal::convert_using_rfunction(case_specific[var_case_specific], "as.matrix");
  Eigen::Map<Eigen::MatrixXd> M_case_specific = as<Eigen::Map<Eigen::MatrixXd> >(case_specific_num);

  Eigen::MatrixXd Matrix_trans((N_cats-1)*N_ind,(N_cats-1)*var_case_specific.length());



  for(int indi = 1 ; indi <= N_ind ; indi++)
  {

    for(int var = 0 ; var <= var_case_specific.length()-1; var++){
      if (is_ref_alt_spe(var) == 1){ // CASE WHERE hinc[car] que tienen ref car
        Eigen::MatrixXd Block_ind =  M_case_specific.block((indi-1) * N_cats, var, 1, 1);
        Eigen::MatrixXd Block_RES = Eigen::kroneckerProduct(Block_ind, ONES_Q1).eval();
        Matrix_trans.block((indi-1) * (N_cats-1), (N_cats-1)*(var), N_cats-1, (N_cats-1)) = Block_RES;
      }else{
        Eigen::MatrixXd Block_ind =  M_case_specific.block((indi-1) * N_cats, var, 1, 1);
        Eigen::MatrixXd Block_RES = Eigen::kroneckerProduct(Block_ind, Iden_Q1).eval();
        Matrix_trans.block((indi-1) * (N_cats-1), (N_cats-1)*(var), N_cats-1, (N_cats-1)) = Block_RES;

      }
    }
  }

  // Create vector relation between category and order in dataset
  CharacterVector ordered_cat = case_specific["alternatives"];
  String cat_1 = ordered_cat[0];
  NumericVector cat_index = NumericVector::create(Named(cat_1 , 0));

  for(int j_1 = 1; j_1 < N_cats ; j_1++){
    cat_1 = ordered_cat[j_1];
    cat_index.push_back(j_1, cat_1);
  }
  int count_re_var = 0;
  for(int ind_b_var = 0; ind_b_var < var_case_specific.length() ; ind_b_var++){

    Eigen::MatrixXd Block_cat = Matrix_trans.block(0,ind_b_var*(N_cats-1),Matrix_trans.rows(), N_cats-1);

    // Rcout << "paso1" << std::endl;

    if( effect_specific_for2[ind_b_var] != "" ){
      String cat_loop = effect_specific_for2[ind_b_var];
      int Var_to_keep ;
      // if(ad_or_ref == "reference"){
      if (is_ref_alt_spe(ind_b_var) == 0){
        // EL INDICE DE LAS COLUMNAS CON LAS QUE ME VOY A QUEDAR PARA LA MATRIX EXTENDIDA
        Var_to_keep = cat_index[cat_loop];
      }else{
        int var1 = cat_index[cat_loop];
        Var_to_keep = var1-1;}
      Matrix_trans.block(0,count_re_var,Matrix_trans.rows(), 1) = Block_cat.block(0, Var_to_keep, Matrix_trans.rows(),1);
      count_re_var = count_re_var+1;

      // }else if(ad_or_ref == "adjacent"){
      //   // Rcout << "es adjc" << std::endl;
      //   if (is_ref_alt_spe(ind_b_var) == 0){
      //     Var_to_keep = cat_index[cat_loop];
      //   }else{
      //     // int var1 = cat_index[cat_loop];
      //     Var_to_keep = N_cats;}
      //   // NumericVector Ind_adj(N_cats+2);
      //   // Ind_adj[Var_to_keep+1] = 1;
      //   // Ind_adj[Var_to_keep] = -1;
      //   // Ind_adj.erase(0);
      //   // Ind_adj.erase(N_cats+1);
      //   // Ind_adj.erase(N_cats);
      //   // NumericVector Ind_rep = rep(Ind_adj,N_ind);
      //   NumericVector to_mul(Block_cat.rows());
      //   NumericVector Adj_vec(Block_cat.rows());
      //   // NumericVector to_mul = wrap(Block_cat.rowwise().sum());
      //   if (is_ref_alt_spe(ind_b_var) == 0){
      //     NumericVector Ind_adj(N_cats+2);
      //     Ind_adj[Var_to_keep+1] = 1;
      //     Ind_adj[Var_to_keep] = -1;
      //     Ind_adj.erase(0);
      //     Ind_adj.erase(N_cats+1);
      //     Ind_adj.erase(N_cats);
      //     NumericVector Ind_rep = rep(Ind_adj,N_ind);
      //     to_mul = wrap(Block_cat.rowwise().sum());
      //     Adj_vec = Ind_rep * to_mul;
      //     // Rcout << "vec resultante 1" << std::endl;
      //     // print(Adj_vec);
      //   }else{
      //     NumericVector Ind_adj(N_cats+2);
      //     Ind_adj[N_cats-1] = -1;
      //     // Ind_adj[Var_to_keep] = -1;
      //     Ind_adj.erase(0);
      //     Ind_adj.erase(N_cats+1);
      //     Ind_adj.erase(N_cats);
      //     // print(Ind_adj);
      //     NumericVector Ind_rep = rep(Ind_adj,N_ind);
      //     to_mul = wrap(-Block_cat.block(0, 0, Matrix_trans.rows(),1));
      //     Adj_vec = Ind_rep * to_mul;
      //     // Rcout << "vec resultante 2" << std::endl;
      //     // print(Adj_vec);
      //   }
      // NumericVector Adj_vec = Ind_rep * to_mul;
      // Eigen::Map<Eigen::VectorXd> Adj_vec1(Rcpp::as<Eigen::Map<Eigen::VectorXd>>(Adj_vec));
      // Matrix_trans.block(0,count_re_var,Matrix_trans.rows(), 1) = Adj_vec1;
      // count_re_var = count_re_var+1;
      // print(Adj_vec);
      // }else{stop("Error: Unknown link. Valid options are reference or adjacent");}

    }else{
      Matrix_trans.block(0,count_re_var,Matrix_trans.rows(), N_cats-1) = Block_cat;
      count_re_var = count_re_var+N_cats-1;
    }
  }



  Eigen::MatrixXd Matrix_trans1 = Matrix_trans.block(0,0,Matrix_trans.rows(), count_re_var);

  // print(head1(Matrix_trans1));

  // Aca esta la parte de lo que seria el diseno complete, asi que si quiero un solo intercepto eliminare las
  // J-1 primeras filas y reemplazare por el vector original de 1

  // String conditional = "conditional";
  if (intercept == "conditional"){
    Matrix_trans1 = Matrix_trans1.block(0, N_cats-2, Matrix_trans1.rows(), Matrix_trans1.cols() - (N_cats-2) );
    Matrix_trans1.col(0) = M_case_specific.col(0);
  }

  // print(head1(Matrix_trans1));

  return Matrix_trans1;
}


Eigen::MatrixXd Extend_All_design(DataFrame Final_mat,
                                  DataFrame Var_alt,
                                  CharacterVector var_alt_specific,
                                  int N_ind, int N_cats,
                                  String ref_cat,
                                  String intercept
                                    // , String ad_or_ref
){

  LogicalVector any_alternative_specific = !is_na(var_alt_specific); // TRUE IF ANY

  CharacterVector var_case_specific = Var_Not_In(Final_mat, var_alt_specific);

  Eigen::MatrixXd Ex_case_M = Extend_case_specific(Final_mat, N_cats, N_ind, var_alt_specific, Var_alt , ref_cat, intercept);
  // Eigen::MatrixXd Ex_case_M = Extend_case_specific(Final_mat, N_cats, N_ind, var_alt_specific, Var_alt , ref_cat, ad_or_ref);
  Eigen::MatrixXd Design_Matrix;

  // print(head1(Ex_case_M));

  if (any_alternative_specific[0]) {
    Eigen::MatrixXd Ex_alt_M = Extend_alt_specific(Final_mat, N_cats, N_ind, var_alt_specific);
    // Eigen::MatrixXd Ex_alt_M = Extend_alt_specific(Final_mat, N_cats, N_ind, var_alt_specific, ad_or_ref);
    int rows_t = Ex_case_M.rows();
    Design_Matrix.conservativeResize(rows_t,Ex_case_M.cols()+Ex_alt_M.cols());
    Design_Matrix.block(0,0,rows_t,Ex_case_M.cols()) = Ex_case_M;
    Design_Matrix.block(0,Ex_case_M.cols(),rows_t,Ex_alt_M.cols()) = Ex_alt_M;
  }else{Design_Matrix = Ex_case_M;}

  return Design_Matrix;
}

List Extend_Response(DataFrame Final_mat ){
  Environment base_env("package:base");
  Function my_unique = base_env["unique"];
  Function my_length = base_env["length"];
  Function my_matrix = base_env["matrix"];
  Function my_asnumeric = base_env["as.numeric"];
  Function my_asfactor = base_env["as.factor"];
  SEXP N_cat_1 = my_length(my_unique(Final_mat["alternatives"]));
  int N_cat = Rcpp::as<int>(N_cat_1);
  NumericVector y11 = my_asnumeric(my_asfactor(Final_mat["choice"])) ;
  DataFrame Y_Ext = my_transpose(my_matrix(y11 ,  _["nrow"] = N_cat));
  // Rcout << "Y_Ext" << std::endl;
  // Rcout << Y_Ext << std::endl;

  NumericMatrix Y_Ext1 = internal::convert_using_rfunction(Y_Ext, "as.matrix");
  Eigen::MatrixXd Y_n2 = as<Eigen::Map<Eigen::MatrixXd> >(Y_Ext1-1);
  // Eigen::MatrixXd Y_n2 = as<Eigen::Map<Eigen::MatrixXd> >(Y_Ext1);
  Y_n2.conservativeResize(Y_n2.rows(), Y_n2.cols() - 1);
  // Rcout << "Y_n2" << std::endl;
  // Rcout << Y_n2 << std::endl;

  return List::create(
    Named("N_cat") = N_cat,
    Named("Y_Ext") = Y_n2
  );
}

List cdf::select_data_nested(Formula formula,
                             String individuals,
                             String Alternatives,
                             CharacterVector ref_cat,
                             CharacterVector var_alt_specific,
                             DataFrame input_data,
                             String intercept,
                             String predict
                               //   ,
                               // String ad_or_ref
) {

  Environment stats_env1("package:utils");
  Function head = stats_env1["head"];


  List Formula_l = formula_entry(formula);
  // print(Formula_l);
  SEXP Var_spe_alt1 = Formula_l["Var_alt"];
  String Response = Formula_l["Response"];
  DataFrame Var_spe_alt = Rcpp::as<DataFrame>(Var_spe_alt1);

  CharacterVector var_informatives = CharacterVector::create("Alternatives",
                                                             "Cat_ref_vec",
                                                             "individuals");
  var_informatives[0] = Alternatives;
  var_informatives[2] = individuals;

  // Rcout << var_informatives << std::endl;

  List Final_mat = All_pre_data(Formula_l["formula_model"],
                                input_data,
                                var_informatives,
                                Response,
                                ref_cat,
                                predict);

  SEXP MA1 = Final_mat["data_output"];;
  DataFrame MA11 = Rcpp::as<DataFrame>(MA1);
  DataFrame Final_mat1 = my_AsNumericMatrix(MA11);
  CharacterVector Levels = Final_mat["Levels"];
  DataFrame new1 = Var_spe_alt;

  LogicalVector any_alternative_specific = !is_na(var_alt_specific);
  CharacterVector colnames_final_m = Var_spe_alt.names();

  if (any_alternative_specific[0]) {
    NumericVector x(colnames_final_m.length());
    for(int indi = 0 ; indi < var_alt_specific.length(); indi++){
      String var_1 = var_alt_specific[indi];
      int indi_var = new1.findName(var_1);
      x[indi_var] = indi_var;
      new1.erase(indi_var);
    }
    colnames_final_m = colnames_final_m[x==0]; // Case where there no alternative specific variables
  }
  DataFrame new2 = new1;
  CharacterVector Names_design;
  Environment base_base("package:base");
  Function my_paste = base_base["paste"];



  for(int indi = 0 ; indi < new2.cols(); indi++){
    CharacterVector colnames = new2.names();
    CharacterVector var_1 = new2[indi];
    String var111 = var_1[0];
    // print(var_1[0]);
    if(var_1[0] != ""){
      String var11 = colnames[indi];
      CharacterVector a1 = my_paste(var11, var111,  _["collapse"] = "");
      Names_design.push_back(a1[0]);
    }else{
      // String conditional = "conditional";
      if (intercept == "conditional" && indi == 0){
        String var11 = colnames[indi];
        CharacterVector a1 = my_paste(var11, var111,  _["collapse"] = "");
        Names_design.push_back(a1[0]);
      }else{
        for(int cats = 0 ; cats < Levels.length()-1; cats++){
          String var11 = colnames[indi];
          String var12 = Levels[cats];
          String var1 = my_paste(var11, var12,  _["collapse"] = "");
          Names_design.push_back(var1);
        }
      }
    }
  }

  if(any_alternative_specific[0]){
    for(int i = 0 ; i < var_alt_specific.length(); i++ ){
      Names_design.push_back(var_alt_specific[i]);
    }
  }

  String ref_cat1 = ref_cat[ref_cat.length()-1];
  List Response_L = Extend_Response(Final_mat1);
  Eigen::MatrixXd Response_M = Response_L["Y_Ext"];

  //   print(Var_spe_alt); // All variables including the intercept
  //   print(var_alt_specific); // JUst the specific for alternatives

  Eigen::MatrixXd Design_Matrix = Extend_All_design(Final_mat1,
                                                    Var_spe_alt,
                                                    var_alt_specific,
                                                    Response_M.rows(),
                                                    Response_L["N_cat"],
                                                              ref_cat1,
                                                              intercept
                                                      // , ad_or_ref
  );


  //
  //   print(head(Names_design));
  //   print(head(Design_Matrix));

  return List::create(
    _["Design_Matrix"] = Design_Matrix,
    _["Response_M"] = Response_M,
    _["Names_design"] = Names_design,
    _["categories_order"] = Levels
  );
}

Eigen::VectorXd Logistic::in_open_corner(const Eigen::VectorXd& p) const
{
  Eigen::VectorXd pi = p;
  int J = pi.size() + 1;
  for(int j=0; j<J-1; ++j)
  { pi[j] = std::max(_epsilon_0, std::min(pi[j], 1-_epsilon_1)); }
  double sum = pi.sum();
  if(sum > 1-_epsilon_1)
  {
    for(int j=0; j<J-1; ++j)
    { pi[j] *= (1.-_epsilon_1)/sum;  }
  }
  return pi;
}

Logistic::Logistic(void) {
}
double Logistic::cdf_logit(const double& value) const
{
  boost::math::logistic dist(0., 1.);
  return boost::math::cdf(dist, value);
}

double Logistic::cdf_logit_complement(const double& value) const
{
  boost::math::logistic dist(0., 1.);
  return boost::math::cdf(complement(dist, value));
}

double Logistic::pdf_logit(const double& value) const
{
  boost::math::logistic dist(0., 1.);
  return boost::math::pdf(dist, value);
}
double Logistic::qdf_logit(const double& value) const
{
  boost::math::logistic dist(0., 1.);
  return boost::math::quantile(dist, value);
}

Eigen::VectorXd Logistic::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  boost::math::logistic dist(0., 1.);
  for (int i = 0; i<=vector.size()-1; i++)
    vector(i) = quantile(dist, vector(i));
  return vector;
}


Normal::Normal(void) {
}
double Normal::cdf_normal(const double& value) const
{
  boost::math::normal norm;
  return boost::math::cdf(norm, value);
}
double Normal::cdf_normal_complement(const double& value) const
{
  boost::math::normal norm;
  return boost::math::cdf(complement(norm, value));
}
double Normal::pdf_normal(const double& value) const
{
  boost::math::normal norm;
  return boost::math::pdf(norm, value);
}
double Normal::qdf_normal(const double& value) const
{
  boost::math::normal norm;
  return boost::math::quantile(norm, value);
}
Eigen::VectorXd Normal::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  boost::math::normal norm;
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = quantile(norm, vector(i));
  return vector;
}

Cauchy::Cauchy(void) {
}

double Cauchy::cdf_cauchy(const double& value) const
{
  double _location = 0.0;
  double _scale = 1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return boost::math::cdf(extreme_value, value);
}
double Cauchy::cdf_cauchy_complement(const double& value) const
{
  double _location = 0.0;
  double _scale = 1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return boost::math::cdf(complement(extreme_value, value));
}
double Cauchy::pdf_cauchy(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return pdf(extreme_value, value);
}
double Cauchy::qdf_cauchy(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  return quantile(extreme_value, value);
}
Eigen::VectorXd Cauchy::InverseLinkQuantileFunction(Eigen::VectorXd vector ){
  double _location = 0.0;
  double _scale =1.0;
  boost::math::cauchy_distribution<> extreme_value(_location, _scale);
  for (int i = 0; i<=vector.rows()-1; i++)
    vector(i) = quantile(extreme_value, vector(i));
  return vector;
}

Student::Student(void) {
}

double Student::cdf_student(const double& value, const double& freedom_degrees) const
{
  double z;
  if(freedom_degrees < 2 * pow(value, 2) )
  { z = boost::math::ibeta(freedom_degrees * 0.5, 0.5, freedom_degrees / (freedom_degrees + pow(value, 2))) * 0.5; }
  else
  { z = boost::math::ibetac(0.5, freedom_degrees * 0.5, pow(value, 2) / (freedom_degrees + pow(value, 2))) * 0.5; }
  if(value > 0)
  { return 1-z; }
  else
  {return z; }
}

// double Student::cdf_student(const double& value, const double& freedom_degrees) const
// {
//   boost::math::students_t_distribution<> student(freedom_degrees);
//   return boost::math::cdf(student,value);
// }
//
// double Student::pdf_student(const double& value, const double& freedom_degrees) const
// {
//   boost::math::students_t_distribution<> student(freedom_degrees);
//   return boost::math::pdf(student,value);
// }

double Student::qdf_student(const double& value, const double& freedom_degrees) const
{
  boost::math::students_t_distribution<> student(freedom_degrees);
  return boost::math::quantile(student,value);
}

double Student::pdf_student(const double& value, const double& freedom_degrees) const
{ return pow( freedom_degrees/(freedom_degrees + pow(value, 2)) ,
              (1+freedom_degrees) * 0.5 ) / ( pow(freedom_degrees,0.5) *
                boost::math::beta(freedom_degrees*0.5, 0.5) ); }


Gumbel::Gumbel(void) {
  // Rcout << "Gumbel is being created" << endl;
}
double Gumbel::cdf_gumbel(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return boost::math::cdf(extreme_value, value);
}
double Gumbel::cdf_gumbel_complement(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return boost::math::cdf(complement(extreme_value, value));
}
double Gumbel::pdf_gumbel(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return pdf(extreme_value, value);
}
double Gumbel::qdf_gumbel(const double& value) const
{
  double _location = 0.0;
  double _scale =1.0;
  boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
  return quantile(extreme_value, value);
}

Gompertz::Gompertz(void) {
  // Rcout << "Gompertz is being created" << endl;
}

double Gompertz::pdf_gompertz(const double& value) const
{
  // double _mu = 0.0;
  // double _sigma = 1.0;

  // return exp(1+value-exp(value)) ; }
  // return exp(value-(exp(value)-1)) ; }
  // return exp(value)*( exp(-(exp(value)-1)) ) ; }
  return (exp(value-exp(value))) ; }

// double Gompertz::cdf_gompertz(const double& value) const
// {
//   double _location = 0.0;
//   double _scale =1.0;
//   boost::math::extreme_value_distribution<> extreme_value(_location, _scale);
//   return boost::math::cdf(complement(extreme_value, -value));
// }
double Gompertz::cdf_gompertz(const double& value) const
{
  // double _mu = 0.0;
  // double _sigma = 1.0;
  // return  1 - exp(-(exp(value)-1)); }
  return(1-exp(-exp(value)));}

double Gompertz::qdf_gompertz(const double& value) const
{
  // double _mu = 0.0;
  // double _sigma = 1.0;
  return  log( -log(1-value)); }

Laplace::Laplace(void) {
  // Rcout << "Laplace is being created" << endl;
}

double Laplace::pdf_laplace(const double& value) const
{
  boost::math::laplace dist(0., 1.);
  return boost::math::pdf(dist, value);
}

double Laplace::cdf_laplace(const double& value) const
{ boost::math::laplace dist(0., 1.);
  return boost::math::cdf(dist, value);
}
double Laplace::cdf_laplace_complement(const double& value) const
{ boost::math::laplace dist(0., 1.);
  return boost::math::cdf(complement(dist, value));
}

double Laplace::qdf_laplace(const double& value) const
{ boost::math::laplace dist(0., 1.);
  return boost::math::quantile(dist, value);
}

Noncentralt::Noncentralt(void) {
  // Rcout << "Laplace is being created" << endl;
}

double Noncentralt::pdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const
{
  boost::math::non_central_t dist(freedom_degrees, non_centrality);
  return boost::math::pdf(dist, value);
}

double Noncentralt::cdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const
{ boost::math::non_central_t dist(freedom_degrees, non_centrality);
  return boost::math::cdf(dist, value);
}
double Noncentralt::cdf_non_central_t_complement(const double& value, const double& freedom_degrees, const double& non_centrality) const
{ boost::math::non_central_t dist(freedom_degrees, non_centrality);
  return boost::math::cdf(complement(dist, value));
}

double Noncentralt::qdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const
{ boost::math::non_central_t dist(freedom_degrees, non_centrality);
  return boost::math::quantile(dist, value);
}


// RCPP_MODULE(exportmod){
//   using namespace Rcpp ;
//   class_<cdf>("cdf")
//     .constructor()
//   ;
// }

// RCPP_MODULE(exportmoddev){
//   using namespace Rcpp ;
//   class_<cdf>("cdf")
//     .constructor()
//   ;
//   class_<Logistic>("Logistic")
//     .derives<cdf>("cdf")
//     .constructor()
//     .method( "InverseLinkCumulativeFunction", &Logistic::InverseLinkCumulativeFunction )
//   ;
// }

