#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class KalmanFilter {
public:

  // state vector
  VectorXd x_;

  // state covariance matrix
  MatrixXd P_;

  // state transistion matrix
  MatrixXd F_;

  // process covariance matrix
  MatrixXd Q_;

  // measurement matrix
  MatrixXd H_;

  // measurement covariance matrix
  MatrixXd R_;

  KalmanFilter();
  virtual ~KalmanFilter();

  void Predict();
  void PrepareForLaserUpdate(const Eigen::VectorXd &z);
  void PrepareForRadarUpdate(const Eigen::VectorXd &z);
  void Update();

 private:
  VectorXd z_;
  VectorXd z_pred_;
};

#endif /* KALMAN_FILTER_H_ */
