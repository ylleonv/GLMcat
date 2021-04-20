#ifndef REFERENCEF_H_
#define REFERENCEF_H_
#include "cdf.h"

class ReferenceF : public virtual Logistic, Normal, Cauchit, Student, Gumbel, Gompertz, Laplace, Noncentralt {
public:

  ReferenceF();

  virtual Eigen::VectorXd inverse_logistic(const Eigen::VectorXd& eta1) const ;
  virtual Eigen::MatrixXd inverse_derivative_logistic(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_normal(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_normal(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_cauchit(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_cauchit(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gompertz(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gompertz(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_gumbel(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_gumbel(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_laplace(const Eigen::VectorXd& eta) const;
  virtual Eigen::MatrixXd inverse_derivative_laplace(const Eigen::VectorXd& eta) const ;

  virtual Eigen::VectorXd inverse_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const;
  // virtual Eigen::MatrixXd student_D(const Eigen::VectorXd& eta, const double& freedom_degrees) const;
  virtual Eigen::MatrixXd inverse_derivative_student(const Eigen::VectorXd& eta, const double& freedom_degrees) const ;

  virtual Eigen::VectorXd inverse_noncentralt(const Eigen::VectorXd& eta, const double& freedom_degrees, const double& mu) const;
  virtual Eigen::MatrixXd inverse_derivative_noncentralt(const Eigen::VectorXd& eta, const double& freedom_degrees, const double& mu) const ;

  // List Discrete_CM(std::string response_categories, std::string individual_choice,
  //                std::string individuals,
  //                StringVector explanatory_global,
  //                StringVector category_specific,
  //                std::string cdf,
  //                StringVector initial_order,
  //                std::string reference_category,
  //                DataFrame dataframe,
  //                std::string design,
  //                double freedom_degrees);

};

#endif
