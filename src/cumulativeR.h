#ifndef CUMULATIVER_H_
#define CUMULATIVER_H_
#include "cdf.h"

class CumulativeR : virtual public cdf, Logistic, Normal, Cauchit, Student, Gumbel, Gompertz, Laplace, Noncentralt{
public:
  CumulativeR();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const;

  virtual Eigen::VectorXd inverse_gompertz(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gompertz(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const;
  virtual Eigen::MatrixXd inverse_derivative_student(const Eigen::VectorXd& eta,const double& freedom_degrees) const ;

  virtual Eigen::VectorXd inverse_laplace(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_laplace(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gumbel(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gumbel(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_noncentralt(const Eigen::VectorXd& eta, const double& freedom_degrees, const double& mu) const;
  virtual Eigen::MatrixXd inverse_derivative_noncentralt(const Eigen::VectorXd& eta, const double& freedom_degrees, const double& mu) const ;


  // List GLMcum(std::string response,
  //             StringVector explanatory_complete,
  //             StringVector explanatory_parallel,
  //             std::string cdf,
  //             SEXP categories_order,
  //             DataFrame dataframe,
  //             StringVector beta_t,
  //             Eigen::VectorXd beta_init);


};

#endif


