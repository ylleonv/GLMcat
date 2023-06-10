#ifndef cdf_H
#define cdf_H

#include <RcppEigen.h>
using namespace std;
using namespace Rcpp;

class cdf{
public:
  double _epsilon_0 = 1e-10;
  double _epsilon_1 = 1e-6;

  std::string concatenate(std::string x, std::string level);

  List All_pre_data_or(Formula formula,
                       DataFrame input_data,
                       CharacterVector categories_order,
                       CharacterVector parallel_effect,
                       std::string threshold = "NA",
                       std::string ratio = "non_cum");

  List All_pre_data_NEWDATA(Formula formula,
                            DataFrame NEWDATA,
                            CharacterVector categories_order,
                            CharacterVector parallel_effect,
                            int N_cats
  );

  List select_data_nested(Formula formula,
                          String individuals,
                          String Alternatives,
                          CharacterVector ref_cat,
                          CharacterVector var_alt_specific,
                          DataFrame input_data,
                          String intercept,
                          String predict
                            //   ,
                            // String ratio
  );

  cdf();
};

class Logistic : virtual public cdf{
public:
  virtual Eigen::VectorXd in_open_corner(const Eigen::VectorXd& p) const;
  virtual double cdf_logit(const double& value) const;
  virtual double cdf_logit_complement(const double& value) const;
  virtual double pdf_logit(const double& value) const;
  virtual double qdf_logit(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Logistic();
};

class Normal : virtual public cdf{
public:
  virtual double cdf_normal(const double& value) const;
  virtual double cdf_normal_complement(const double& value) const;
  virtual double pdf_normal(const double& value) const;
  virtual double qdf_normal(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Normal();
};

class Cauchy : virtual public cdf{
public:
  virtual double cdf_cauchy(const double& value) const;
  virtual double cdf_cauchy_complement(const double& value) const;
  virtual double pdf_cauchy(const double& value) const;
  virtual double qdf_cauchy(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Cauchy();
};

class Student :  virtual public cdf{
public:
  virtual double cdf_student(const double& value, const double& freedom_degrees) const;
  // virtual double cdf_student_complement(const double& value, const double& freedom_degrees) const;
  virtual double pdf_student(const double& value, const double& freedom_degrees) const;
  virtual double qdf_student(const double& value, const double& freedom_degrees) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Student();
};


class Noncentralt :  virtual public cdf{
public:
  virtual double cdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const;
  virtual double cdf_non_central_t_complement(const double& value, const double& freedom_degrees, const double& non_centrality) const;
  virtual double pdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const;
  virtual double qdf_non_central_t(const double& value, const double& freedom_degrees, const double& non_centrality) const;

  Noncentralt();
};

class Gumbel :  virtual public cdf{
public:
  virtual double cdf_gumbel(const double& value) const;
  virtual double cdf_gumbel_complement(const double& value) const;
  virtual double pdf_gumbel(const double& value) const;
  virtual double qdf_gumbel(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Gumbel();
};

class Gompertz : virtual public cdf{
public:
  virtual double cdf_gompertz(const double& value) const;
  // virtual double cdf_gompertz_complement(const double& value) const;
  virtual double pdf_gompertz(const double& value) const;
  virtual double qdf_gompertz(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Gompertz();
};

class Laplace : virtual public cdf{
public:
  virtual double cdf_laplace(const double& value) const;
  virtual double cdf_laplace_complement(const double& value) const;
  virtual double pdf_laplace(const double& value) const;
  virtual double qdf_laplace(const double& value) const;

  Eigen::VectorXd InverseLinkQuantileFunction(Eigen::VectorXd vectordis);

  Laplace();
};


#endif
