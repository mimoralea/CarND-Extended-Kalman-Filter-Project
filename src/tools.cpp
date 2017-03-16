#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}
Tools::~Tools() {}

/**
 * Calculates the error between the
 * estimated values after the Kalman cycle
 * and the ground truth we had
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
        VectorXd rmse(4);
        rmse << 0,0,0,0;

        // make sure we have the same number of measurements
        if(estimations.size() != ground_truth.size()
           || estimations.size() == 0) {
                cout << "Invalid estimation or ground_truth data" << endl;
                return rmse;
        }

        // calculate the squares differences
        for(unsigned int i=0; i < estimations.size(); ++i) {
                VectorXd residual = estimations[i] - ground_truth[i];
                residual = residual.array() * residual.array();
                rmse += residual;
        }

        // average and take square root
        rmse = rmse / estimations.size();
        rmse = rmse.array().sqrt();
        return rmse;
}

/**
 * Calculate Jacobian Matrix used to
 * linearize the measurements from RADAR
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

        // extract state values
        float px = x_state(0), py = x_state(1), vx = x_state(2), vy = x_state(3);
        float px2_p_py2 = px * px + py * py;

        // protect too small values from passing
        if(fabs(px2_p_py2) < 0.0001) {
                throw 1;
        }

        // calculate the values of the derivatives
        MatrixXd Hj(3,4);
        Hj(0, 0) = px / sqrt(px2_p_py2);
        Hj(0, 1) = py / sqrt(px2_p_py2);
        Hj(0, 2) = 0;
        Hj(0, 3) = 0;

        Hj(1, 0) = -py / px2_p_py2;
        Hj(1, 1) = px / px2_p_py2;
        Hj(1, 2) = 0;
        Hj(1, 3) = 0;

        Hj(2, 0) = py * (vx * py - vy * px) / pow(px2_p_py2, 3/2.);
        Hj(2, 1) = px * (vy * px - vx * py) / pow(px2_p_py2, 3/2.);
        Hj(2, 2) = px / sqrt(px2_p_py2);
        Hj(2, 3) = py / sqrt(px2_p_py2);

        // return Jacobian matrix
        return Hj;
}
