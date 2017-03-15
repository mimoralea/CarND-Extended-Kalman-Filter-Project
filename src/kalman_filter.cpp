#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update() {
  VectorXd y = z_ - z_pred_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::PrepareForLaserUpdate(const VectorXd &z) {
  z_pred_ = H_ * x_;
  z_ = z;
}

void KalmanFilter::PrepareForRadarUpdate(const VectorXd &z) {
  float px = x_[0], py = x_[1], vx = x_[2], vy = x_[3];

  float rho = sqrt(px*px + py*py);
  float phi = atan(py/px);
  float rho_dot = (px*vx + py*vy)/sqrt(px*px + py*py);

  z_pred_ = VectorXd(3);
  z_pred_ << rho, phi, rho_dot;
  z_ = z;
}
