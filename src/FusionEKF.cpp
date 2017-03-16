#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Constructor initializes matrices with default values
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

        // noise for the process covariance
        sigma_ax = 9, sigma_ay = 9;

        // state matrix
        ekf_.x_ = VectorXd(4);
        // state covariance matrix P
        ekf_.P_ = MatrixXd(4, 4);

        //the initial transition matrix F_
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;
}

/**
 * Destructor
 */
FusionEKF::~FusionEKF() {}

/**
 * Kalman filter initial state and processing
 */
void FusionEKF::Init(const MeasurementPackage &measurement_pack) {

        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
        if (is_initialized_) {
                // only initialize once
                return;
        }

        float px, py, vx = 0, vy = 0;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

                // extract the RADAR measurements and convert from
                // Polar to Cartesian coordinates
                float range = measurement_pack.raw_measurements_[0];
                float bearing = measurement_pack.raw_measurements_[1];
                float range_rate = measurement_pack.raw_measurements_[2];

                // calculate position and velocity
                px = range * cos(bearing);
                py = range * sin(bearing);
                vx = range_rate * cos(bearing);
                vy = range_rate * sin(bearing);
                // note that we are using the velocity measurements on this
                // because radar can measure velocity (unlike laser)

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

                // if it is laser, just grab the raw x, y coordinates
                px = measurement_pack.raw_measurements_[0];
                py = measurement_pack.raw_measurements_[1];
        }

        // set the state and state covariance matrices
        ekf_.x_ << px, py, vx, vy;
        ekf_.P_ << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;

        // ensure we mark the timestamp
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
}

/**
 * Mark timestamp and calculate elapsed
 * seconds since last measurement
 */
float FusionEKF::ProcessTimestamp(const long timestamp) {
        float dt = (timestamp - previous_timestamp_) / 1000000.0;
        previous_timestamp_ = timestamp;
        return dt;
}

/**
 * Add the elapsed time since last
 * measurement into the state transition matrix
 */
void FusionEKF::ProcessStateTransMatrix(const float dt) {
        ekf_.F_(0, 2) = dt;
        ekf_.F_(1, 3) = dt;
}

/**
 * Calculate the Qv and G matrix to get the
 * process covariance matrix `Q`
 */
void FusionEKF::ProcessCovarianceMatrix(const float dt) {
        MatrixXd Qv = MatrixXd(2, 2);
        Qv << sigma_ax, 0,
                0, sigma_ay;

        MatrixXd G = MatrixXd(4, 2);
        G << dt * dt / 2, 0,
                0, dt * dt / 2,
                dt, 0,
                0, dt;

        // calculate process covariance matrix
        ekf_.Q_ = G * Qv * G.transpose();
}

/**
 * Processes the measurement packages as they enter the
 * Kalman Filter implementation.
 * It initializes, Predicts and Updates measurements
 * from both RADAR and LASER data effectively fusing them
 */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

        /*****************************************************************************
         *  Initialization
         ****************************************************************************/
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
        VectorXd z_pred;
        switch (measurement_pack.sensor_type_) {
        case MeasurementPackage::RADAR:
                try {
                        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
                        ekf_.R_ = R_radar_;
                        z_pred = ekf_.PrepareForRadarUpdate();
                } catch(...) {
                        // Jacobian error with a very small value -- ignore current update
                        return;
                }
                break;
        case MeasurementPackage::LASER:
                ekf_.H_ = H_laser_;
                ekf_.R_ = R_laser_;
                z_pred = ekf_.PrepareForLaserUpdate();
                break;
        default:
                break;
        }

        // update either RADAR or LASER measurement
        ekf_.Update(measurement_pack.raw_measurements_,
                    z_pred);
}
