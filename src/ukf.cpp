#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::endl;

// used to compute the RMSE later
Tools tools;
/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .6;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;
  
  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  P_ << 9., 0, 0, 0, 0,
        0, 9., 0, 0, 0,
        0, 0, .9, 0, 0,
        0, 0, 0, 3., 0,
        0, 0, 0, 0, 3.;
  
  limit_ = 2 * n_aug_ + 1;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;
  delta_t_;
  Xsig_pred_ = MatrixXd(n_x_, limit_);
  weights_ = VectorXd(limit_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if( !use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
    return;
  
  if( !use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
    return;
  
  // cout << "ProcessMeasurement" << endl;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {   
    time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phy = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      
      float x = rho * cos(phy);
      float y = rho * sin(phy);
      x_ << x, y, 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    //set weights
    double wt_n = lambda_ + n_aug_;
    weights_(0) = lambda_/wt_n;
    for(int i=1; i < limit_; i++)
    {
        weights_(i) = 0.5/(wt_n);
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
   
  /*********************************************************
                            PREDICTION STEP
  *********************************************************/
  delta_t_ = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_= meas_package.timestamp_;
  
  Prediction();
  
  /*********************************************************
                            UPDATE STEP
  *********************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  
  // cout << x_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t_ the change in time (in seconds) between the last
 * measurement and this one.
 */
 
void UKF::Prediction() {
  // cout << "Prediction" << endl;
  MatrixXd Xsig_aug = GenerateSigmaPoints();
  
  PredictSigmaPoints(Xsig_aug);
    
  PredictMeanCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // cout << "UpdateLidar" << endl;
  VectorXd z_pred;
  MatrixXd S;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig;
  PredictLaserMeasurement(&z_pred, &S, &Zsig);
  // cout << "after predictlaser meas method" << endl;
  //create example vector for incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < limit_; i++)
  {
	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  VectorXd z_diff = Zsig.col(i) - z_pred;
	  
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  // cout << "Ater TC calculation" << endl;
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  
  x_ = x_ + K * z_diff;
  
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // cout << "UpdateRadar" << endl;
  VectorXd z_pred;
  MatrixXd S;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig;
  PredictRadarMeasurement(&z_pred, &S, &Zsig);
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_);
  z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i = 0; i < limit_; i++)
  {
	  VectorXd x_diff = Xsig_pred_.col(i) - x_;
	  VectorXd z_diff = Zsig.col(i) - z_pred;
	  
	  //normalize angles
    
	  x_diff(3) = tools.NormalizeAngle(x_diff(3));
	  z_diff(1) = tools.NormalizeAngle(z_diff(1));
	  
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  z_diff(1) = tools.NormalizeAngle(z_diff(1));
  
  x_ = x_ + K * z_diff;
  
  P_ = P_ - K * S * K.transpose();
}

MatrixXd UKF::GenerateSigmaPoints()
{
  // cout << "GenerateSigmaPoints" << endl;
  /*********************************************************
          Generate sigma points
  *********************************************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, limit_);
  
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++){
	  MatrixXd cal = sqrt(lambda_ + n_aug_) * L.col(i);
	  Xsig_aug.col(i + 1) = x_aug + cal;
	  Xsig_aug.col(i + 1 + n_aug_) = x_aug - cal;
  }
  
  return Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_aug)
{
  // cout << "PredictSigmaPoints" << endl;
  /*********************************************************
          Predict sigma points
  *********************************************************/
  //predict sigma points
  for (int i = 0; i < limit_; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);
	
	//predicted state values
    double px_p, py_p;
	  double v_p = v;
    double yaw_p = yaw + yawd*delta_t_;
    double yawd_p = yawd;
	
    if (fabs(yawd) > 0.001) //non zero
    {
      px_p = p_x + (v/yawd * (sin(yaw + yawd * delta_t_) - sin(yaw) ) );
      py_p = p_y + (v/yawd * (-cos(yaw + yawd*delta_t_) + cos(yaw) ) );
    }
    else
    {
      px_p = p_x + v*cos(yaw)*delta_t_;
      py_p = p_y + v*sin(yaw)*delta_t_;
    }
    
    //add noise
    double delta_t_2by2 = 0.5*delta_t_*delta_t_;
    px_p += delta_t_2by2*cos(yaw)*nu_a;
    py_p += delta_t_2by2*sin(yaw)*nu_a;
    v_p += delta_t_ * nu_a;
    yaw_p += delta_t_2by2*nu_yawdd;
    yawd_p += delta_t_*nu_yawdd;
    
    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
   }  
}

void UKF::PredictMeanCovariance()
{  
  // cout << "PredictMeanCovariance" << endl;
  //predict state mean
  x_.fill(0.0);
  for(int i=0; i<limit_; i++)
  { //iterate over sigma points
      x_ += weights_(i) * Xsig_pred_.col(i);
  }
  P_.fill(0.0);
  //predict state covariance matrix
  for(int i=0; i<limit_; i++)
  {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      //Normalize radian angle
      x_diff(3) = tools.NormalizeAngle(x_diff(3));
      
      P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  // cout << "PredictRadarMeasurement" << endl;
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  MatrixXd Zsig = MatrixXd(n_z_, limit_);
  
  S.fill(0.0);
  
  z_pred.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < limit_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                          //rho
    Zsig(1,i) = atan2(p_y,p_x);                                   //phi
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v ) / Zsig(0,i);   //rho_dot
    
	//calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
 
  for (int i = 0; i < limit_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools.NormalizeAngle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_ * std_radr_ , 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0,std_radrd_ * std_radrd_;
  S += R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  // cout << "PredictRadarMeasurement" << endl;
  //mean predicted measurement
  int n_z = 2;
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd Zsig = MatrixXd(n_z, limit_);
  
  S.fill(0.0);
  
  z_pred.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < limit_; i++) {  //2n+1 simga points
    Zsig(0,i)  = Xsig_pred_(0,i);
    Zsig(1,i)  = Xsig_pred_(1,i);

	//calculate mean predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }
  
  //calculate innovation covariance matrix S
 
  for (int i = 0; i < limit_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << std_laspx_ * std_laspx_ , 0,
       0, std_laspy_ * std_laspy_;
  S += R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}