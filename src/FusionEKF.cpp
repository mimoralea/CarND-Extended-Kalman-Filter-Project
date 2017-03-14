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
  R_radar_ = MatrixXd(3, 3);

  H_laser_ = MatrixXd(2, 4);

  Hj_ = MatrixXd(3, 4);

	//create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;

	//measurement covariance
	ekf_.R_ = MatrixXd(2, 2);
	ekf_.R_ << 0.0225, 0,
             0, 0.0225;

	//measurement matrix
	ekf_.H_ = MatrixXd(2, 4);
	ekf_.H_ << 1, 0, 0, 0,
             0, 1, 0, 0;

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

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "radar" << endl;
      float ro = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float ro_dot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << ro * cos(phi), ro * sin(phi), ro_dot * cos(phi), ro_dot * sin(phi);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "lidar" << endl;
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "initialized" << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

	//compute the time elapsed between the current and previous measurements
  float sigma_ax = 9;
  float sigma_ay = 9;
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  cout << dt << endl;

  //1. Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

  cout << ekf_.F_ << endl;

	//2. Set the process covariance matrix Q
	MatrixXd Qv = MatrixXd(2, 2);
	Qv << sigma_ax, 0,
        0, sigma_ay;

  MatrixXd G = MatrixXd(4, 2);
  G << dt * dt / 2, 0,
       0, dt * dt / 2,
       dt, 0,
       0, dt;
  ekf_.Q_ = G * Qv * G.transpose();
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
