
#include "ukf.h"

#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Main constructor, initlaizing variable
 */
UKF::UKF() {

 /* default state of teh filter*/
 is_initialized_ = false;

 /*Time stamp in microsecond*/
// time_us_ = 0;
 //state dimension for radar
 n_x_ = 5;

 //augmented state dimension for ladar to support unscented operation
 n_aug_ = 7;
 // if this is false, laser measurements will be ignored (except during init)
 use_laser_ = true;

 // if this is false, radar measurements will be ignored (except during init)
 use_radar_ = true;

 // initial state vector
 x_ = VectorXd(n_x_);

 // initial covariance matrix
 P_ = MatrixXd(n_x_, n_x_);

 // Process noise standard deviation (longitudinal acceleration in m/s^2)
 std_a_ = 2.15; //original: 30 is too high when compared with other vals

 // Process noise standard deviation (yaw acceleration in rad/s^2)
 std_yawdd_ = 2; //original: 30 is too too high

 // Laser measurement noise standard deviation position1 in m
 std_laspx_ = 0.15;

 // Laser measurement noise standard deviation position2 in m
 std_laspy_ = 0.15;

 // Radar measurement noise standard deviation radius in m
 std_radr_ = 0.3;

 // Radar measurement noise standard deviation angle in rad
 std_radphi_ = 0.03;

 // Radar measurement noise standard deviation radius change in m/s
 std_radrd_ = 0.3;



 //lambda spreading parameter for unscented operations
 lambda_ = 3 - n_x_;

 //Corresponding weights for UKF
 weights_ = VectorXd(  (2 * n_aug_) + 1);

 /* Initial weight value is different */
 weights_(0) = lambda_  / ( lambda_ + n_aug_ );

 /* Rest of the weights from formula in notes*/
 for (int i = 1; i<2 * n_aug_ + 1; i++)
	 weights_(i) = 1 / (2 * (lambda_ + n_aug_));


 // Prediction covariance matrix, can be initialised in constructor
  P_ <<   1, 0, 0, 0, 0,
		  0, 1, 0, 0, 0,
		  0, 0, 1, 0, 0,
		  0, 0, 0, 1, 0,
		  0, 0, 0, 0, 1;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

 /* If very first packet, then initialise */

 if (!is_initialized_) {

	 // Initialize the state ekf_.x_
	 x_ << 1, 1, 1, 1, 0;


	 if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

	 // Use the first measurement to initialise UKF

	 double range = meas_package.raw_measurements_[0];   //magnitude

	 double phi = tools.NormAngle(meas_package.raw_measurements_[1]);   //angle in radians, start a check on the initial value too

	 double range_dot = meas_package.raw_measurements_[2];   //range rate

	 double px = range * cos(phi);
	 double py = range * sin(phi);
	 double vx = range_dot * cos(phi);
	 double vy = range_dot * sin(phi);

	 double val = atan2(vy, vx);

	 double phi_dot = tools.NormAngle(val);

	 //if initial values are zero, start from non-zero values for faster convergence
	 if (px == 0 || py == 0) {
		 px =   0.001;
		 py =   0.001;
	 }
	 x_ << px, py, 1, 1, 1;

	 }
	 else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
	 // Set state x_ to the first measurement.
		 x_(0) = meas_package.raw_measurements_[0];
		 x_(1) = meas_package.raw_measurements_[1];

		 //if initial values are zero
		 /*if (px == 0 && py == 0)
		 {
		 px = py = 0.001;
		 }
*/
		 //x_ << px, py, 0, 0, 0;
	 }

	 /* save the first timestam as it is*/
	 timestamp_prev_ = meas_package.timestamp_;

	 // done initializing
	 is_initialized_ = true;
	 return;


   } //End of first initialization


 /* If not first time then do normal KF  Prediction  and Update steps   */


  //compute the time elapsed between the current and previous measurements
 double delta_t = (meas_package.timestamp_ - timestamp_prev_) / 1000000.0;	//dt - expressed in seconds
 timestamp_prev_ = meas_package.timestamp_;

 // make prediction of sigma points
 Prediction(delta_t);


 //update radar measurememt and state estimation
 if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
	 UpdateRadar(meas_package);

 //update laser measurememt and state estimation
 else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
	 UpdateLidar(meas_package);


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
 /**
 Estimate the object's location. Modify the state
 vector, x_. Predict sigma points, the state, and the state covariance matrix.
 */

 //First make required matrices for augmented space


 VectorXd x_aug = VectorXd(n_aug_);


 MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
 P_aug.setZero();   // initialize to 0

 //sigma point matrix, 7 * 15
 MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


 //augmented state
 x_aug.head(5) = x_;   // take first 5 values from previous a posterioir state, rest are  0
 x_aug(5) = 0;
 x_aug(6) = 0;

 //augmented P matrix
 P_aug.fill(0.0);
 P_aug.topLeftCorner<5, 5>() = P_;   //fill the top left section
 P_aug(5,5) = std_a_*std_a_;
 P_aug(6,6) = std_yawdd_*std_yawdd_;
 //P_aug.bottomRightCorner<2, 2>() = Q;

 //State Matrix
 MatrixXd A = P_aug.llt().matrixL();

 //augmented sigma points
 //set first column of sigma point matrix

 Xsig_aug.col(0) = x_aug;

 //set remaining sigma points
 for (int i = 0; i < n_aug_; i++){
	 Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
	 Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
 }


 //Augmented Sigma point prediction
 Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

 //predict sigma points
 for (int i = 0; i< 2 * n_aug_ + 1; i++){
	 //Obtaining individual Xaug elements
	 double p_x = Xsig_aug(0, i);
	 double p_y = Xsig_aug(1, i);
	 double v = Xsig_aug(2, i);
	 double yaw = Xsig_aug(3, i);
	 double yawd = Xsig_aug(4, i);
	 double nu_a = Xsig_aug(5, i);
	 double nu_yawdd = Xsig_aug(6, i);


	 double px_p, py_p;  //predicted state values

	 //avoid division by zero, code from lectures
	 if (fabs(yawd) > 0.001)
	 	 {
		 px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
		 py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
	 	 }
	 else
	 	 {
		px_p = p_x + v*delta_t*cos(yaw);
		py_p = p_y + v*delta_t*sin(yaw);
	 	 }

	 double v_p = v;
	 double yaw_p = yaw + yawd*delta_t;
	 double yawd_p = yawd;

	 //adding  noise related parameters
	 px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
	 py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
	 v_p = v_p + nu_a*delta_t;
	 yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
	 yawd_p = yawd_p + nu_yawdd*delta_t;

	 //write predicted sigma point into right column
	 Xsig_pred_(0, i) = px_p;
	 Xsig_pred_(1, i) = py_p;
	 Xsig_pred_(2, i) = v_p;
	 Xsig_pred_(3, i) = yaw_p;
	 Xsig_pred_(4, i) = yawd_p;
    } //end of the long for loop


 //predict state mean
 x_.fill(0.0);


 //iteratation with sigma points
 for (int i = 0; i < 2 * n_aug_ + 1; i++)
	 x_ = x_ + weights_(i) * Xsig_pred_.col(i);


 //predicted state covariance matrix
 P_.fill(0.0);
 //iterate over sigma points
 for (int i = 0; i < 2 * n_aug_ + 1; i++)
    {

	 VectorXd x_diff = Xsig_pred_.col(i) - x_;  	 // state difference
	 x_diff(3) = tools.NormAngle(x_diff(3));    //Normalize angle before use
	 P_ = P_ + weights_(i) * x_diff * x_diff.transpose();  //set the rest of Pmatrix with weights
   }


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
 /**
 Use lidar data to update the belief about the object's
 position. Modify the state vector, x_, and covariance, P_.
 Calculate the lidar NIS.
 */

 // sensor state dimension
 int n_z = 2;

 //create matrix for sigma points in measurement space
 MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

 //transform sigma points into measurement space
 for (int i = 0; i < 2 * n_aug_ + 1; i++){

	 double p_x = Xsig_pred_(0, i);
	 double p_y = Xsig_pred_(1, i);
	 double v = Xsig_pred_(2, i);
	 double yaw = Xsig_pred_(3, i);

	 // measurement model
	 Zsig(0, i) = Xsig_pred_(0, i); //px
	 Zsig(1, i) = Xsig_pred_(1, i); //py

 }

 //mean predicted measurement
 VectorXd z_pred = VectorXd(n_z);
 z_pred.fill(0.0);

 for (int i = 0; i < 2 * n_aug_ + 1; i++)
	 z_pred = z_pred + weights_(i) * Zsig.col(i);


 //measurement covariance matrix S
 MatrixXd S = MatrixXd(n_z, n_z);
 S.fill(0.0);
 for (int i = 0; i < 2 * n_aug_ + 1; i++) {

	 VectorXd z_diff = Zsig.col(i) - z_pred;

	 z_diff(1) = tools.NormAngle(z_diff(1)); //  normalization of angle before use

	 S = S + weights_(i) * z_diff * z_diff.transpose();
 }

 //add measurement noise covariance matrix
 MatrixXd R = MatrixXd(n_z, n_z);
 R << std_laspx_*std_laspx_, 0,
 0, std_laspy_*std_laspy_;
 S = S + R;


 // cross correlation G
 MatrixXd G = MatrixXd(n_x_, n_z);
 G.fill(0.0);

 //calculate cross correlation matrix
 for (int i = 0; i < 2 * n_aug_ + 1; i++) {

	 VectorXd z_diff = Zsig.col(i) - z_pred;

	 z_diff(1) = tools.NormAngle(z_diff(1));   //Norm angle before use

	 // state difference
	 VectorXd x_diff = Xsig_pred_.col(i) - x_;
	 //Norm angle before use
	 x_diff(3) = tools.NormAngle(x_diff(3));

	 //Cross-correlation matrix
	 G = G + weights_(i) * x_diff * z_diff.transpose();
 }

 //Kalman gain K;
 MatrixXd K = G * S.inverse();

 //actual measurement
 VectorXd z = VectorXd(n_z);
 z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];// 0.0, 0.0;

 //residual
 VectorXd z_diff = z - z_pred;

 //Norm angle before use
 z_diff(1) = tools.NormAngle(z_diff(1));

 //Final kalman update
 x_ = x_ + K * z_diff;
 P_ = P_ - K *    S * K.transpose();

 //Calculate the lidar NIS.
 nis_lidar_= z_diff.transpose() * S.inverse() * z_diff;


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
 /**
 Use radar data to update the belief about the object's
 position. Modify the state vector, x_, and covariance, P_.
 Calculate the radar NIS.
 */

 // sensor state dimension
 int n_z = 3;

 //create matrix for sigma points in measurement space
 MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

 //transform sigma points into measurement space
 for (int i = 0; i < 2 * n_aug_ + 1; i++)  { //2n+1 simga points


	 double p_x = Xsig_pred_(0, i);
	 double p_y = Xsig_pred_(1, i);
	 double v = Xsig_pred_(2, i);
	 double yaw = Xsig_pred_(3, i);

	 double v1 = cos(yaw)*v;
	 double v2 = sin(yaw)*v;

	 // measurement model
	 Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y); //range
	 Zsig(1, i) = atan2(p_y, p_x);         //phi
	 Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); //range_dot
 }

 //mean predicted measurement
 VectorXd z_pred = VectorXd(n_z);
 z_pred.fill(0.0);
 for (int i = 0; i < 2 * n_aug_ + 1; i++)
	 z_pred = z_pred + weights_(i) * Zsig.col(i);


 //measurement covariance matrix S
 MatrixXd S = MatrixXd(n_z, n_z);
 S.fill(0.0);
 for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

	 VectorXd z_diff = Zsig.col(i) - z_pred;   //difference

	 //Normalize angle before use
	 z_diff(1) = tools.NormAngle(z_diff(1));

	 S = S + weights_(i) * z_diff * z_diff.transpose();
 }

 //Obtain R matrix
 MatrixXd R = MatrixXd(n_z, n_z);
 R << std_radr_*std_radr_, 0, 0,
 0, std_radphi_*std_radphi_, 0,
 0, 0, std_radrd_*std_radrd_;
 S = S + R;


 // cross correlation matrix G
 MatrixXd G = MatrixXd(n_x_, n_z);
 G.fill(0.0);

 //calculate cross correlation matrix
 for (int i = 0; i < 2 * n_aug_ + 1; i++) { //2n+1 simga points

	 //residual
	 VectorXd z_diff = Zsig.col(i) - z_pred;
	 //Norm angle before use
	 z_diff(1) = tools.NormAngle(z_diff(1));

	 // state difference
	 VectorXd x_diff = Xsig_pred_.col(i) - x_;
	 //Norm angle before use
	 x_diff(3) = tools.NormAngle(x_diff(3));

	 G = G + weights_(i) * x_diff * z_diff.transpose();
 }

 //Gain matrix , Kalman Filter
 MatrixXd K = G * S.inverse();

 //actual measurement
 VectorXd z = VectorXd(n_z);
 z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

 //error
 VectorXd z_diff = z - z_pred;

 //Norm angle before use
 z_diff(1) = tools.NormAngle(z_diff(1));


 x_ = x_ + K * z_diff;   // kalman update
 P_ = P_ - K * S * K.transpose();

 nis_radar_ = z_diff.transpose() * S.inverse() * z_diff;

}
