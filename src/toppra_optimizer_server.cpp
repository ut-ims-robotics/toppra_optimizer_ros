#include "ros/ros.h"

#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/parametrizer/const_accel.hpp>
#include <toppra/toppra.hpp>

#include <toppra/solver/qpOASES-wrapper.hpp>

#include <toppra/solver/glpk-wrapper.hpp>

#include <trajectory_msgs/JointTrajectory.h>
#include <trajectory_msgs/JointTrajectoryPoint.h>

#include <toppra_optimizer_ros/OptimizeAction.h>
#include <actionlib/server/simple_action_server.h>

#include <std_msgs/Header.h>

#include <numeric>
#include <vector>
#include <string>

void transpose(std::vector<std::vector<double> > &b)
{
  if (b.size() == 0)
    return;

  std::vector<std::vector<double> > trans_vec(b[0].size(), std::vector<double>());

  for (int i = 0; i < b.size(); i++)
  {
    for (int j = 0; j < b[i].size(); j++)
    {
      trans_vec[j].push_back(b[i][j]);
    }
  }

  b = trans_vec;    // <--- reassign here
}

void formatVecToMat(const std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> &vec, Eigen::MatrixXd &mat)
{
  mat.resize(vec.at(0).rows(), vec.size());
  for (size_t i = 0; i < vec.size(); i++)
  {
    mat.col(i) = vec.at(i);
  }
}


double length(const std::vector<double> &q0, const std::vector<double> &q1)
{
  auto dof = q0.size();
  double retval = 0;

  for (size_t i = 0; i < dof; i++) {
    retval += (q0[i] - q1[i]) * (q0[i] - q1[i]);
  }
  return sqrt(retval) / dof;
}

inline toppra::Vectors makeVectors(
    const std::vector<std::vector<toppra::value_type>> &v)
{
  toppra::Vectors ret;
  for (auto vi : v) {
    toppra::Vector vi_eigen(vi.size());
    for (std::size_t i = 0; i < vi.size(); i++) {
      vi_eigen(i) = vi[i];
    }
    ret.push_back(vi_eigen);
  }
  return ret;
}


inline toppra::BoundaryCond makeBoundaryCond(
    int order, const std::vector<toppra::value_type> &values)
{
  toppra::BoundaryCond cond;
  cond.order = order;
  cond.values.resize(values.size());
  for (std::size_t i = 0; i < values.size(); i++) {
    cond.values(i) = values[i];
  }
  return cond;
}

toppra_optimizer_ros::OptimizeFeedback feedback_;
toppra_optimizer_ros::OptimizeResult result_;

typedef actionlib::SimpleActionServer<toppra_optimizer_ros::OptimizeAction> Server;

void execute(const toppra_optimizer_ros::OptimizeGoalConstPtr& goal, Server* as)
{

  trajectory_msgs::JointTrajectory trajectory_to_optimize = goal -> trajectory;
  std::cout << trajectory_to_optimize << std::endl;
  
  int plan_size = trajectory_to_optimize.points.size();

  std::vector<std::vector<double>> joint_positions;
  for (int i = 0; i < plan_size; i++)
  {
    joint_positions.push_back(trajectory_to_optimize.points[i].positions);
  }

  // create piecewise polynomial geometric path (taken from the file in the issue. line 419)
  toppra::Vector times;
  int dof = trajectory_to_optimize.joint_names.size();
  int numPoints = joint_positions.size();
  times.resize(numPoints);
  times(0) = 0.;
  for (int i = 1; i < numPoints; i++)
  {
    times(i) = times(i - 1) + length(joint_positions[i - 1], joint_positions[i]);
  }

  //boundary condition. 2 means we specify acceleration, after that starting/ending acceleration is specified
  toppra::BoundaryCond bc = makeBoundaryCond(2, { 0, 0, 0, 0, 0, 0, 0 });
  std::array<toppra::BoundaryCond, 2> bc_type{ bc, bc };

  auto toppra_joint_positions = makeVectors(joint_positions);

  //combine into path
  auto path = std::make_shared<toppra::PiecewisePolyPath>(toppra_joint_positions, times, bc_type);

  std::vector<int> points(joint_positions.size());
  std::iota (std::begin(points), std::end(points), 0); 

  //set vel and acc limit and put them as constrains

  //////////////////////////////////////////////////
  //
  // IMPORTANT!!!
  //
  // WHAT DO WE DO WITH THE VELOCITY AND ACCELERATION LIMITS?
  //
  //////////////////////////////////////////////////
  auto velLimitLower = -0.21750 * toppra::Vector::Ones(dof);
  auto velLimitUpper = 0.21750 * toppra::Vector::Ones(dof);
  auto accLimitLower = -0.18750 * toppra::Vector::Ones(dof);
  auto accLimitUpper = 0.18750 * toppra::Vector::Ones(dof);

  toppra::LinearConstraintPtr ljv, lja;
  ljv = std::make_shared<toppra::constraint::LinearJointVelocity>
    (velLimitLower, velLimitUpper);
  lja = std::make_shared<toppra::constraint::LinearJointAcceleration>
    (accLimitLower, accLimitUpper);


  toppra::LinearConstraintPtrs constrains =
    toppra::LinearConstraintPtrs{ljv, lja};

  //make an algorithm with path and constrains
  toppra::PathParametrizationAlgorithmPtr algo =
    std::make_shared<toppra::algorithm::TOPPRA>(constrains, path);

  //compute the algorithm and get the resulted data
  toppra::ReturnCode rc1 = algo->computePathParametrization(0, 0);
  toppra::ParametrizationData pd = algo->getParameterizationData();

  toppra::Vector gridpoints =
      pd.gridpoints; // Grid-points used for solving the discretized problem.
  toppra::Vector vsquared =
      pd.parametrization; // Output parametrization (squared path velocity)

  std::shared_ptr<toppra::parametrizer::ConstAccel> ca =
      std::make_shared<toppra::parametrizer::ConstAccel>(path, gridpoints,
                                                          vsquared);

  //time(?) length of the path
  Eigen::Matrix<toppra::value_type, 1, 2> interval2;
  interval2 = ca->pathInterval();

  // if (true)
  //     std::cout << "ca->validate() = " << ca->validate() << std::endl;

  //how many thetas you want
  const int length2 = plan_size; //same amount as was in the original plan
  toppra::Vector times2 =
      toppra::Vector::LinSpaced(length2, interval2(0), interval2(1)); //make timestamps to get corresponding thetas

    std::cout << "interval2 = " << interval2 << std::endl;

  //take thetas, acc and vel for computed timestamps
  toppra::Vectors path_pos2;
  path_pos2 = ca->eval(times2, 0);
  toppra::Vectors path_vel2;
  path_vel2 = ca->eval(times2, 1);
  toppra::Vectors path_acc2;
  path_acc2 = ca->eval(times2, 2);

  //converting the data to matrix format
  Eigen::MatrixXd path_pos2_ = Eigen::MatrixXd::Zero(dof, length2);
  Eigen::MatrixXd path_vel2_ = Eigen::MatrixXd::Zero(dof, length2);
  Eigen::MatrixXd path_acc2_ = Eigen::MatrixXd::Zero(dof, length2);
  formatVecToMat(path_pos2, path_pos2_);
  formatVecToMat(path_vel2, path_vel2_);
  formatVecToMat(path_acc2, path_acc2_);

  //rotate the matrix to be the shape of [plan_size][dof]
  Eigen::MatrixXd path_pos2_proper = Eigen::MatrixXd::Zero(length2, dof);
  Eigen::MatrixXd path_vel2_proper = Eigen::MatrixXd::Zero(length2, dof);
  Eigen::MatrixXd path_acc2_proper = Eigen::MatrixXd::Zero(length2, dof);
  path_pos2_proper = path_pos2_.reverse().transpose().reverse();
  path_vel2_proper = path_vel2_.reverse().transpose().reverse();
  path_acc2_proper = path_acc2_.reverse().transpose().reverse();
  times2 = times2.transpose();

  //create plan and start to fill it with computed values
  trajectory_msgs::JointTrajectory final_plan;

  for (int i = 0; i < length2; i++)
  {
    trajectory_msgs::JointTrajectoryPoint plan_increment;
    for (int j = 0; j < dof; j++)
    {
      plan_increment.accelerations.push_back(path_acc2_proper(i,j));
      plan_increment.positions.push_back(path_pos2_proper(i,j));
      plan_increment.velocities.push_back(path_vel2_proper(i,j));
    }
    plan_increment.time_from_start = ros::Duration(times2(i));
    final_plan.points.push_back(plan_increment);
  }

  final_plan.joint_names = trajectory_to_optimize.joint_names;
  final_plan.header = trajectory_to_optimize.header;

  feedback_.optimized_trajectory = final_plan;
  result_.optimization_success = true;

  as->setSucceeded();

}


int main(int argc, char **argv)
{
  ros::init(argc, argv, "toppra_optimizer_action_client");
  ros::NodeHandle _nh;

  Server server(_nh, "optimize_trajectory", boost::bind(&execute, _1, &server), false);
  
  server.start();
  ros::spin();

  return 0;
}

