#include <iostream>
#include "tools.h"
#include <stdexcept>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
    
  if(estimations.size() == 0)
  {
      throw std::logic_error("empty estimations.");
  }
  if(estimations.size() != ground_truth.size() )
  {
      throw std::logic_error("Estimations vector length should be equal to ground truth vectors.");
  }
	//accumulate squared residuals
  VectorXd residual;
  for(int i=0; i < estimations.size(); ++i){
        residual = ground_truth[i] - estimations[i];
        residual = residual.array()*residual.array();
        rmse += residual;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}


float Tools::NormalizeAngle(float angle){
	return atan2(sin(angle), cos(angle) );
}