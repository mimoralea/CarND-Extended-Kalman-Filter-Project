#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}
KalmanFilter::~KalmanFilter() {}

/**
 * Predicted process covariance matrix
 * predicts the new process covariance
 * matrix using the values from the
 * previous iteration
 */
void KalmanFilter::Predict() {
        x_ = F_ * x_;
        MatrixXd Ft = F_.transpose();
        P_ = F_ * P_ * Ft + Q_;
}

/**
 * Calculate Kalman Gain Matrix to be
 * used for adjusting the measurement that
 * we get from the positions and velocities
 */
MatrixXd KalmanFilter::CalculateKalmanGain() {
        MatrixXd Ht = H_.transpose();
        MatrixXd S = H_ * P_ * Ht + R_;
        MatrixXd Si = S.inverse();
        MatrixXd PHt = P_ * Ht;
        MatrixXd K = PHt * Si;
        return K;
}

/**
 * Core update function that works regardless
 * of the sensor type that is being
 * currently evaluated
 */
void KalmanFilter::Update(const VectorXd &z, const VectorXd &z_pred) {

        // calculate Kalman Gain
        MatrixXd K = CalculateKalmanGain();

        // predicted error
        VectorXd y = z - z_pred;

        // calculate new state matrix
        // using Kalman gain to adjust
        // for uncertainty
        x_ = x_ + (K * y);
        long x_size = x_.size();
        MatrixXd I = MatrixXd::Identity(x_size, x_size);

        // Finally update the state covariance matrix
        P_ = (I - K * H_) * P_;
}

/**
 * Adjust the measurements to process a LASER update
 */
VectorXd KalmanFilter::PrepareForLaserUpdate() {

        // calculate predicted value
        // by only using the H matrix
        VectorXd z_pred = H_ * x_;
        return z_pred;
}

/**
 * Adjust the measurements to process a RADAR update
 */
VectorXd KalmanFilter::PrepareForRadarUpdate() {

        // extract the raw values from the state
        float px = x_[0], py = x_[1], vx = x_[2], vy = x_[3];

        // convert Cartesian back to Polar
        float range = sqrt(px*px + py*py);
        float bearing = atan(py/px);
        float range_rate = (px*vx + py*vy)/sqrt(px*px + py*py);

        // predicted state vector
        VectorXd z_pred = VectorXd(3);
        z_pred << range, bearing, range_rate;
        return z_pred;
}
