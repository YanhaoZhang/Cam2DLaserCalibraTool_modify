//
// Created by heyijia on 18-12-11.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>

#include "ceres/ceres.h"
#include "LaseCamCalCeres.h"
#include "pose_local_parameterization.h"
#include "utilities.h"

class PointInPlaneFactor: public ceres::SizedCostFunction<1,7>
{
private:
    Eigen::Vector4d planar_;
    Eigen::Vector3d point_;

public:
    PointInPlaneFactor(Eigen::Vector4d planar,    // planar
    Eigen::Vector3d point)   // point
    : planar_(planar),point_(point)
    {}

    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
};

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d &q)
{
    Eigen::Matrix3d ans;
    ans << 0.0, -q(2), q(1),
            q(2), 0.0, -q(0),
            -q(1), q(0), 0.0;
    return ans;
}
bool PointInPlaneFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {
    Eigen::Vector3d tcl(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond qcl(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d pt_c = qcl.toRotationMatrix() * point_ + tcl;
    residuals[0] = planar_.head(3).transpose() * pt_c + planar_[3];

    if (jacobians)
    {
        if (jacobians[0]) {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);
            Eigen::Matrix<double, 1, 6> jaco_i;
            jaco_i.leftCols<3>() = planar_.head(3);
            jaco_i.rightCols<3>() = planar_.head(3).transpose() * (-qcl.toRotationMatrix() * skewSymmetric(point_));

            jacobian_pose_i.leftCols<6>() = jaco_i;
            jacobian_pose_i.rightCols<1>().setZero();
        }


    }

    return true;
}

void CalibrationTool_SavePlanePoints(const std::vector<Oberserve> obs, const Eigen::Matrix4d Tcl, const std::string path)
{

    std::ofstream fs_planar(path + "planar.txt");
    std::ofstream fs_points(path + "RoiPoints.txt");
    std::ofstream fs_lines(path + "RoiPtOnLines.txt");
    fs_planar << std::setprecision(3);
    fs_points << std::setprecision(3);
    fs_lines << std::setprecision(3);

    for(size_t i = 0; i< obs.size(); ++i) {
        Oberserve obi = obs[i];

        /// save planar
        // transform planar in tag frame to cam frame
        //https://stackoverflow.com/questions/7685495/transforming-a-3d-plane-using-a-4x4-matrix
        Eigen::Vector4d planar_tag(0, 0, 1, 0);  // tag 坐标系下的平面方程
        Eigen::Matrix4d Tctag = Eigen::Matrix4d::Identity();
        Tctag.block(0, 0, 3, 3) = obi.tagPose_Qca.toRotationMatrix();
        Tctag.block(0, 3, 3, 1) = obi.tagPose_tca;
        Eigen::Vector4d planar_cam = (Tctag.inverse()).transpose() * planar_tag;
        fs_planar<< i <<" "<< planar_cam[0] <<" "<< planar_cam[1] <<" "<< planar_cam[2]<<" "<< planar_cam[3]<< std::endl;

        /// save points
        for (size_t j = 0; j < obi.points.size(); ++j) {
            Eigen::Vector4d pt;
            pt << obi.points[j],1;
            pt = Tcl * pt;

            fs_points<< i <<" " << pt[0] <<" "<< pt[1] <<" "<< pt[2]<< std::endl;
        }

        /// save lines
        for (size_t k = 0; k < obi.points_on_line.size(); ++k) {
            Eigen::Vector4d pt;
            pt << obi.points_on_line[k],1;
            pt = Tcl * pt;

            fs_lines<< i <<" " << pt[0] <<" "<< pt[1] <<" "<< pt[2]<< std::endl;
        }
    }

}

/**
 * @brief 利用激光点落在平面上这个约束构建误差方程
 * @param obs 观察到的激光数据
 * @param Tcl transform from laser frame to camera frame
 * @param use_linefitting_data 是否只用两个激光点落在平面上作为约束，一帧激光很多激光点落在平面上，但实际约束相当于只有两个点。
 * @param use_boundary_constraint 是否利用激光点落在标定板两侧边界上这个约束，这个约束非常有效，但是需要使用者知道标定板尺寸和大小。
 *
 */
#define LOSSFUNCTION
void CamLaserCalibration(const std::vector<Oberserve> obs, Eigen::Matrix4d &Tcl, bool use_linefitting_data, bool use_boundary_constraint)
{
    Eigen::Quaterniond q(Tcl.block<3,3>(0,0));

    ceres::Problem problem;
    Eigen::VectorXd pose(7);
    pose << Tcl(0,3),Tcl(1,3),Tcl(2,3),q.x(),q.y(),q.z(), q.w();   // laser to camera (origin)


    for(size_t i = 0; i< obs.size(); ++i)
    {
        Oberserve obi = obs[i];
        // transform planar in tag frame to cam frame
//        https://stackoverflow.com/questions/7685495/transforming-a-3d-plane-using-a-4x4-matrix
        Eigen::Vector4d planar_tag(0,0,1,0);  // tag 坐标系下的平面方程
        Eigen::Matrix4d Tctag = Eigen::Matrix4d::Identity();
        Tctag.block(0,0,3,3) = obi.tagPose_Qca.toRotationMatrix();
        Tctag.block(0,3,3,1) = obi.tagPose_tca;
        Eigen::Vector4d planar_cam = (Tctag.inverse()).transpose() * planar_tag;

        std::vector<Eigen::Vector3d> calibra_pts;
        if(use_linefitting_data)
            calibra_pts =  obi.points_on_line;
        else
            calibra_pts =  obi.points;

        for (size_t j = 0; j < calibra_pts.size(); ++j) {
            Eigen::Vector3d pt =  calibra_pts[j];

            PointInPlaneFactor *costfunction = new PointInPlaneFactor(planar_cam, pt);

#ifdef LOSSFUNCTION
            //ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
            ceres::LossFunction * loss_function = new ceres::CauchyLoss(0.1);
            problem.AddResidualBlock(costfunction, loss_function, pose.data());
#else
            problem.AddResidualBlock(costfunction, NULL, pose.data());
#endif
        }

        // boundary plane constraint
        // 边界约束：每帧激光只利用两个激光点时才用边界约束，不然边界约束权重太小。
        if(use_boundary_constraint && use_linefitting_data)
        {
            // plane pi from ith obs in ith camera frame, 计算标定板边界和相机光心所构建的平面
            // 标定板边界三个定点在标定板坐标系中的坐标

//            length of check board
//            x1 195.2
//            b  135.6
//            a  452.0
//            x2 146.0
//
//            y1 266.8
//            y2 152.0

            Eigen::Vector3d orig(0.01952+0.01356,0.02668+0.01356,0);
            Eigen::Vector3d p1m( 0, 0, 0);
            Eigen::Vector3d p2m( 0.40024, 0, 0);
            Eigen::Vector3d p3m( 0., 0.408, 0);
            p1m -= orig;
            p2m -= orig;
            p3m -= orig;
            // 旋转到相机坐标系
            Eigen::Vector3d p1c = obi.tagPose_Qca.toRotationMatrix() * p1m +  obi.tagPose_tca;
            Eigen::Vector3d p2c = obi.tagPose_Qca.toRotationMatrix() * p2m +  obi.tagPose_tca;
            Eigen::Vector3d p3c = obi.tagPose_Qca.toRotationMatrix() * p3m +  obi.tagPose_tca;

            // 得到平面方程
            Eigen::Vector4d pi1 = pi_from_ppp(p1c, p2c,Eigen::Vector3d( 0, 0, 0 ));
            Eigen::Vector4d pi2 = pi_from_ppp(p1c, p3c,Eigen::Vector3d( 0, 0, 0 ));

            Eigen::Vector3d pt1 =  obi.points.at(0);
            Eigen::Vector3d pt2 =  obi.points.at(obi.points.size()-1);
//            std::cout << " " <<pt1.dot(pi1.head(3))<<std::endl;
            PointInPlaneFactor *costfunction1 = new PointInPlaneFactor(pi1, pt1);
            PointInPlaneFactor *costfunction2 = new PointInPlaneFactor(pi2, pt2);

#ifdef LOSSFUNCTION
            //ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
            ceres::LossFunction * loss_function = new ceres::CauchyLoss(0.1);
            problem.AddResidualBlock(costfunction1, loss_function, pose.data());
            problem.AddResidualBlock(costfunction2, loss_function, pose.data());

#else
            problem.AddResidualBlock(costfunction1, NULL, pose.data());
            problem.AddResidualBlock(costfunction2, NULL, pose.data());
#endif
        }
    }

//    ceres::LocalParameterization* quaternionParameterization = new ceres::QuaternionParameterization;
//    problem.SetParameterization(q_coeffs,quaternionParameterization);
    ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
    problem.AddParameterBlock(pose.data(), 7, local_parameterization);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 100;

    ceres::Solver::Summary summary;
    ceres::Solve (options, &problem, & summary);

    std::cout << summary.FullReport() << std::endl;

    q = Eigen::Quaterniond(pose[6],pose[3],pose[4],pose[5]);

    Tcl.block<3,3>(0,0) = q.toRotationMatrix();
    Tcl.block<3,1>(0,3) << pose[0],pose[1],pose[2];

    /// =============================  analysis code ==============================
    /// Get Information matrix from ceres, used to analysis the Gauge of the system
    Eigen::MatrixXd H(6,6);
    Eigen::MatrixXd b(6,1);
    H.setZero();
    b.setZero();
    double chi = 0;
    for(size_t i = 0; i< obs.size(); ++i)
    {
        Oberserve obi = obs[i];
        // transform planar in tag frame to cam frame
//        https://stackoverflow.com/questions/7685495/transforming-a-3d-plane-using-a-4x4-matrix
        Eigen::Vector4d planar_tag(0,0,1,0);  // tag 坐标系下的平面方程
        Eigen::Matrix4d Tctag = Eigen::Matrix4d::Identity();
        Tctag.block(0,0,3,3) = obi.tagPose_Qca.toRotationMatrix();
        Tctag.block(0,3,3,1) = obi.tagPose_tca;
        Eigen::Vector4d planar_cam = (Tctag.inverse()).transpose() * planar_tag;

        std::vector<Eigen::Vector3d> calibra_pts;
        if(use_linefitting_data)
            calibra_pts =  obi.points_on_line;
        else
            calibra_pts =  obi.points;

        for (size_t j = 0; j < calibra_pts.size(); ++j) {
            Eigen::Vector3d pt =  calibra_pts[j];

            double *res = new double[1];
            double **jaco = new double *[1];
            jaco[0] = new double[1 * 7];

            PointInPlaneFactor *costfunction = new PointInPlaneFactor(planar_cam, pt);

            costfunction->Evaluate(std::vector<double *>{pose.data()}.data(), res, jaco);
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose_i(jaco[0]);
            Eigen::Map<Eigen::Matrix<double, 1, 1>> resd(res);

//            std::cout << jacobian_pose_i << std::endl;
            H += jacobian_pose_i.leftCols<6>().transpose() * jacobian_pose_i.leftCols<6>();
            b -= jacobian_pose_i.leftCols<6>().transpose() * resd;

            chi += resd * resd;

//#ifdef LOSSFUNCTION
//            //ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
//            ceres::LossFunction * loss_function = new ceres::CauchyLoss(0.5);
//            problem.AddResidualBlock(costfunction, loss_function, q_coeffs, t_coeffs);
//#else
//            problem.AddResidualBlock(costfunction, NULL, q_coeffs, t_coeffs, &scale);
//#endif


        }
    }

//    std::cout << H << std::endl;
    std::cout << "----- H singular values--------:\n";
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    std::cout << svd.singularValues() <<std::endl;
    std::cout <<"\nrecover chi2: " <<chi / 2. << std::endl;

}
/////////////////////////////////////////////////////////////////
struct LineFittingResidfual {
    LineFittingResidfual(double x, double y)
            : x_(x), y_(y) {}

    template <typename T>
    bool operator()(const T* const m, T* residual) const {
        residual[0] = m[0] * T(x_) + m[1] * T(y_) + T(1.);
        return true;
    }

private:
    // Observations for a sample.
    const double x_;
    const double y_;
};

void LineFittingCeres(const std::vector<Eigen::Vector3d> Points, Eigen::Vector2d & Line)
{
    double line[3] = {Line(0),Line(1)};

    ceres::Problem problem;
    for(size_t i = 0; i< Points.size(); ++i)
    {
        Eigen::Vector3d obi = Points[i];

        ceres::CostFunction * costfunction =
                new ceres::AutoDiffCostFunction<LineFittingResidfual,1,2>(
                        new LineFittingResidfual(obi.x(),obi.y()));

#ifdef LOSSFUNCTION
            //ceres::LossFunctionWrapper* loss_function(new ceres::HuberLoss(1.0), ceres::TAKE_OWNERSHIP);
            ceres::LossFunction * loss_function = new ceres::CauchyLoss(0.05);
            problem.AddResidualBlock(costfunction, loss_function, line);
#else
            problem.AddResidualBlock(costfunction, NULL, line);
#endif
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 10;

    ceres::Solver::Summary summary;
    ceres::Solve (options, &problem, & summary);

    Line(0) = line[0];
    Line(1) = line[1];

}