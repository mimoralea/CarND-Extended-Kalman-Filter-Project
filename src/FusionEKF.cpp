#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  // H_radar is the Jacobian or x_

  sigma_ax = 9, sigma_ay = 9;

	//create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Init(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (is_initialized_) {
    return;
  }

  float px, py, vx = 0, vy = 0;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rho_dot = measurement_pack.raw_measurements_[2];

    px = rho * cos(phi);
    py = rho * sin(phi);
    vx = rho_dot * cos(phi);
    vy = rho_dot * sin(phi);

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    px = measurement_pack.raw_measurements_[0];
    py = measurement_pack.raw_measurements_[1];
  }

  ekf_.x_ << px, py, vx, vy;

  previous_timestamp_ = measurement_pack.timestamp_;
  is_initialized_ = true;
}

float FusionEKF::ProcessTimestamp(const long timestamp) {
	float dt = (timestamp - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = timestamp;
  return dt;
}

void FusionEKF::ProcessStateTransMatrix(const float dt) {
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;
}

void FusionEKF::ProcessCovarianceMatrix(const float dt) {
	MatrixXd Qv = MatrixXd(2, 2);
	Qv << sigma_ax, 0,
    0, sigma_ay;

  MatrixXd G = MatrixXd(4, 2);
  G << dt * dt / 2, 0,
    0, dt * dt / 2,
    dt, 0,
    0, dt;
  ekf_.Q_ = G * Qv * G.transpose();
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // initialization
  FusionEKF::Init(measurement_pack);

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  float dt = ProcessTimestamp(measurement_pack.timestamp_);
  ProcessStateTransMatrix(dt);
  ProcessCovarianceMatrix(dt);
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  switch (measurement_pack.sensor_type_) {
  case MeasurementPackage::RADAR:
    try {
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.PrepareForRadarUpdate(measurement_pack.raw_measurements_);
    } catch(...) {
      // Jacobian error - ignore current update
      return;
    }
    break;
  case MeasurementPackage::LASER:
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.PrepareForLaserUpdate(measurement_pack.raw_measurements_);
    break;
  default:
    break;
  }

  ekf_.Update();
}
