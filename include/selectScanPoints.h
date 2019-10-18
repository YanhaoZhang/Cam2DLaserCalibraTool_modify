//
// Created by heyijia on 18-12-11.
//

#ifndef PROJECT_SELECTSCANPOINTS_H
#define PROJECT_SELECTSCANPOINTS_H

#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>

class GetRect
{
public:
    std::vector<cv::Rect> rects;
    bool bDraw;
    cv::Rect r;
    cv::Point base;

    cv::Mat layer;
    cv::Mat working;

public :
    void CallBackFunc(int event, int x, int y);
    static void onMouse(int event, int x, int y, int flags, void* userdata);
    void gettingROI(cv::Mat img);
};

std::vector< Eigen::Vector3d > AutoGetLinePts(const std::vector<Eigen::Vector3d> points, float angle_increment, bool debug = true);
std::vector<cv::Rect> PointsToImg(const std::vector<Eigen::Vector3d> points,bool selectPoint = true);
std::vector<Eigen::Vector3d> GetROIScanPoints(const std::vector<Eigen::Vector3d> points, const std::vector<cv::Rect> rects);


//yanhao
void AutoGetLinePts_useReflector(std::vector<Eigen::Vector3d> points_yh, std::vector<float> points_intensity_yh,
                                 std::vector<Eigen::Vector3d> Points,    std::vector<float> Intensity_oneScan,
                                 float camera_fov, float angle_increment, bool debug);

#endif //PROJECT_SELECTSCANPOINTS_H
