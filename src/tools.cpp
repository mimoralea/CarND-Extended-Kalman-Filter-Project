#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

using namespace std;
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  for(unsigned int i=0; i < estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);

  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float px2_p_py2 = px * px + py * py;

  if(fabs(px2_p_py2) < 0.0001) {
    throw 1;
  }

  Hj(0, 0) = px / sqrt(px2_p_py2);
  Hj(0, 1) = py / sqrt(px2_p_py2);

  Hj(1, 0) = -py / px2_p_py2;
  Hj(1, 1) = px / px2_p_py2;

  Hj(2, 0) = py * (vx * py - vy * px) / pow(px2_p_py2, 3/2.);
  Hj(2, 1) = px * (vy * px - vx * py) / pow(px2_p_py2, 3/2.);

  Hj(2, 2) = px / sqrt(px2_p_py2);
  Hj(2, 3) = py / sqrt(px2_p_py2);

  return Hj;
}
