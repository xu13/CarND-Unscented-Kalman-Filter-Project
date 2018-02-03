#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

// using namespace std;
// using Eigen::MatrixXd;
// using Eigen::VectorXd;
// using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = Eigen::VectorXd(5);

  // initial covariance matrix
  P_ = Eigen::MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  // DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  // DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
   * TODO:
   *
   * Complete the initialization. See ukf.h for other member properties.
   *
   * Hint: one or more values initialized above might be wildly off...
   */
  n_x_        = 5;
  n_aug_      = 7;
  n_sig_      = 2 * n_aug_ + 1;
  lambda_     = 3 - n_x_;
  weights_    = Eigen::VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (size_t i = 1; i < n_sig_; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }
  Xsig_pred_ = Eigen::MatrixXd(n_x_, n_sig_);

  // Measurement noise covariance matrix
  R_rad_       = Eigen::MatrixXd::Zero(3, 3);
  R_rad_(0, 0) = std_radr_ * std_radr_;
  R_rad_(1, 1) = std_radphi_ * std_radphi_;
  R_rad_(2, 2) = std_radrd_ * std_radrd_;

  R_las_ = Eigen::MatrixXd::Zero(2, 2);
  R_las_.fill(0.0);
  R_las_(0, 0) = std_laspx_ * std_laspx_;
  R_las_(1, 1) = std_laspy_ * std_laspy_;

  is_initialized_     = false;
  previous_timestamp_ = 0;
}

UKF::~UKF(){ }

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO:
   *
   * Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      x_ << rho * std::cos(phi), rho * std::sin(phi), 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    }
    P_.fill(0.0);
    P_(0, 0) = 1; // sigma_px^2
    P_(1, 1) = 1; // sigma_py^2
    P_(2, 2) = 1; // sigma_v^2
    P_(3, 3) = 1; // sigma_yaw^2
    P_(4, 4) = 1; // sigma_yawd^2

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_     = true;
    return;
  }

  // Prediction
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; // dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  // Update
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
  } else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    // Lidar updates
    UpdateLidar(meas_package);
  }
} // UKF::ProcessMeasurement

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
   * TODO:
   *
   * Complete this function! Estimate the object's location. Modify the state
   * vector, x_. Predict sigma points, the state, and the state covariance matrix.
   */
  // Augment state vector
  Eigen::VectorXd x_aug(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5)      = 0;
  x_aug(6)      = 0;

  // Augment covariance matrix
  Eigen::MatrixXd Q(2, 2);
  Q << std_a_ * std_a_, 0, 0, std_yawdd_ * std_yawdd_;
  Eigen::MatrixXd P_aug(n_aug_, n_aug_);
  P_aug << P_, Eigen::MatrixXd::Zero(n_x_, 2), Eigen::MatrixXd::Zero(2, n_x_), Q;

  // Gerneate sigma points
  Eigen::MatrixXd L = P_aug.llt().matrixL(); // square root of P_aug
  Eigen::MatrixXd Xsig_aug(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;
  for (size_t i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + std::sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - std::sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // Predict sigma points
  for (size_t i = 0; i < n_sig_; i++) {
    const Eigen::VectorXd& x = Xsig_aug.col(i);
    float px       = x(0);
    float py       = x(1);
    float v        = x(2);
    float yaw      = x(3);
    float yawd     = x(4);
    float nu_a     = x(5);
    float nu_yawdd = x(6);

    Eigen::VectorXd f(n_x_);
    Eigen::VectorXd g(n_x_);
    g << 0.5 * delta_t * delta_t * std::cos(yaw) * nu_a,
      0.5 * delta_t * delta_t * std::sin(yaw) * nu_a,
      delta_t * nu_a,
      0.5 * delta_t * delta_t * nu_yawdd,
      delta_t * nu_yawdd;

    if (std::abs(yawd) <= 1e-3) {
      f << v * std::cos(yaw) * delta_t,
        v * std::sin(yaw) * delta_t,
        0,
        yawd * delta_t,
        0;
    } else {
      f << v / yawd * (std::sin(yaw + yawd * delta_t) - std::sin(yaw)),
        v / yawd * (-std::cos(yaw + yawd * delta_t) + std::cos(yaw)),
        0,
        yawd * delta_t,
        0;
    }
    Xsig_pred_.col(i) = x.head(n_x_) + f + g;
  }

  // Predict state mean
  x_.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance matrix
  P_.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    // State difference
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
} // UKF::Prediction

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   * TODO:
   *
   * Complete this function! Use lidar data to update the belief about the object's
   * position. Modify the state vector, x_, and covariance, P_.
   *
   * You'll also need to calculate the lidar NIS.
   */
  // Transform sigma points into measurement space
  unsigned int n_z = 2;
  Eigen::MatrixXd Zsig(n_z, n_sig_);
  for (size_t i = 0; i < n_sig_; i++) {
    const Eigen::VectorXd& x = Xsig_pred_.col(i);
    float px = x(0);
    float py = x(1);
    Zsig.col(i) << px, py;
  }

  // Calculate mean predicted measurement
  Eigen::VectorXd z_pred(n_z);
  z_pred.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  Eigen::MatrixXd S(n_z, n_z);
  S.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    // Residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  S += R_las_;

  //
  // UKF Update
  //
  // Calculate cross correlation matrix
  Eigen::MatrixXd Tc(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) { // 2n+1 simga points
    // Residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

    // State difference
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  Eigen::MatrixXd K = Tc * S.inverse();

  // Residual
  Eigen::VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // Update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // Calculate NIS
  float epsilon = z_diff.transpose() * S.inverse() * z_diff;
  nis_las_.push_back(epsilon);
} // UKF::UpdateLidar

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * TODO:
   *
   * Complete this function! Use radar data to update the belief about the object's
   * position. Modify the state vector, x_, and covariance, P_.
   *
   * You'll also need to calculate the radar NIS.
   */
  // Transform sigma points into measurement space
  unsigned int n_z = 3;
  Eigen::MatrixXd Zsig(n_z, n_sig_);
  for (size_t i = 0; i < n_sig_; i++) {
    const Eigen::VectorXd& x = Xsig_pred_.col(i);
    float px  = x(0);
    float py  = x(1);
    float v   = x(2);
    float yaw = x(3);

    float rho  = std::sqrt(px * px + py * py);
    float phi  = std::atan2(py, px);
    float phid = (px * std::cos(yaw) * v + py * std::sin(yaw) * v) / rho;
    Zsig.col(i) << rho, phi, phid;
  }

  // Calculate mean predicted measurement
  Eigen::VectorXd z_pred(n_z);
  z_pred.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  Eigen::MatrixXd S(n_z, n_z);
  S.fill(0.0);
  for (size_t i = 0; i < n_sig_; i++) {
    // Residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  S += R_rad_;

  //
  // UKF Update
  //
  // Calculate cross correlation matrix
  Eigen::MatrixXd Tc(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) { // 2n+1 simga points
    // Residual
    Eigen::VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    // State difference
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K
  Eigen::MatrixXd K = Tc * S.inverse();

  // Residual
  Eigen::VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

  // Angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2. * M_PI;

  // Update state mean and covariance matrix
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // Calculate NIS
  float epsilon = z_diff.transpose() * S.inverse() * z_diff;
  nis_rad_.push_back(epsilon);
} // UKF::UpdateRadar
