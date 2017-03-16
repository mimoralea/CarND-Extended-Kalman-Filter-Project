#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;

class FusionEKF {
public:
        KalmanFilter ekf_;

        FusionEKF();
        virtual ~FusionEKF();

        void ProcessMeasurement(const MeasurementPackage &measurement_pack);

private:
        bool is_initialized_;
        long previous_timestamp_;
        float sigma_ax, sigma_ay;

        Tools tools;
        MatrixXd R_laser_;
        MatrixXd R_radar_;
        MatrixXd H_laser_;

        void Init(const MeasurementPackage &measurement_pack);

        void ProcessStateTransMatrix(const float dt);
        void ProcessCovarianceMatrix(const float dt);
        float ProcessTimestamp(const long timestamp);
};

#endif /* FusionEKF_H_ */
