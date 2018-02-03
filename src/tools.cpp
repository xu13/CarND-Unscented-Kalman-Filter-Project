#include <iostream>
#include "tools.h"

Tools::Tools(){ }

Tools::~Tools(){ }

Eigen::VectorXd Tools::CalculateRMSE(const std::vector<Eigen::VectorXd>& estimations,
                                     const std::vector<Eigen::VectorXd>& ground_truth)
{
  /**
   * TODO:
   * Calculate the RMSE here.
   */
  size_t N = estimations.size();
  assert(N != 0);
  assert(N == ground_truth.size());

  Eigen::VectorXd rmse(estimations[0].size());
  rmse.fill(0.0);
  for (size_t i = 0; i < N; i++) {
    Eigen::VectorXd e = estimations[i] - ground_truth[i];
    e     = e.array() * e.array();
    rmse += e;
  }
  rmse /= N;
  rmse  = rmse.array().sqrt();

  return rmse;
}
