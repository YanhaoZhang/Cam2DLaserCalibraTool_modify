#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <sensor_msgs/image_encodings.h>
#include <sensor_msgs/LaserScan.h>

#include <rosbag/bag.h>
#include <rosbag/view.h>

#include "../camera_models/include/EquidistantCamera.h"
#include "../camera_models/include/PinholeCamera.h"
#include "config.h"
#include "utilities.h"
#include "selectScanPoints.h"
#include "calcCamPose.h"
#include "LaseCamCalCeres.h"

template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
  std::cout << name <<std::endl;
  T ans;
  if (n.getParam(name, ans))
  {
    ROS_INFO_STREAM("Loaded " << name << ": " << ans);
  }
  else
  {
    ROS_ERROR_STREAM("Failed to load " << name);
    n.shutdown();
  }
  return ans;
}


int main(int argc, char **argv){
  ros::init(argc, argv, "LaserCamCal");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");

  std::string config_file = "/home/yanhao/LaserCameraCal_ws/src/CamLaserCalibraTool/config/calibra_config_pinhole.yaml";
//    std::string config_file;
//    config_file = readParam<std::string>(pnh, "config_file");

    // Yanhao debug
    std::cout<< "Yanhao: Check:----- \n" <<std::endl;
    std::cout<< config_file <<std::endl;




  readParameters(config_file);

  rosbag::Bag bag_input;
  bag_input.open(bag_path, rosbag::bagmode::Read);
  std::vector<std::string> topics;
  topics.push_back(scan_topic_name);
  topics.push_back(img_topic_name);
  rosbag::View views(bag_input, rosbag::TopicQuery(topics));

  // Load apritag pose
  std::vector < CamPose > tagpose;
  LoadCamPoseFromTxt(savePath + "apriltag_pose.txt",tagpose);
  std::cout << "Load apriltag pose size: " << tagpose.size() <<std::endl;

  if(tagpose.size() < 10)
  {
    std::cout << "apriltag pose less than 10." << std::endl;
    return 0;
  }

  /// 移除一些标定板没咋动的数据，防止静止的标定板数据占的比重较多，影响标定结果。比如标定板静止了 1 min, 那会采集非常多的数据，但实际有效的只有一帧数据。
  std::vector< CamPose > sparseTagpose;
  CamPose older = tagpose.at(0);
  sparseTagpose.push_back(older);
  double dist_min = 0.20;  //10cm
  double theta_min = 3.1415926 * 10/ 180.;
  for (int j = 1; j < tagpose.size(); ++j) {
    CamPose newer = tagpose.at(j);
    double dist = (older.twc - newer.twc).norm();
    double theta = 2 * std::acos( ((older.qwc.inverse() * newer.qwc)).w() );
    if( (dist > dist_min) || (fabs(theta) > theta_min) )
    {
      older = newer;
      sparseTagpose.push_back(older);
    }
  }
  tagpose = sparseTagpose;

  // 准备标定数据
  std::vector<Oberserve> obs;
  // 处理激光数据
  int ii_cnt = 10;
  for(rosbag::MessageInstance const m: views){
    if (m.getTopic() == scan_topic_name)
    {
      sensor_msgs::LaserScan::Ptr scan = m.instantiate<sensor_msgs::LaserScan>();
      std::vector<Eigen::Vector3d> Points;
      TranScanToPoints(*scan,Points);



//      ii_cnt++;
//      if(ii_cnt % 20 != 0) continue;

      double timestamp = scan->header.stamp.toSec();
      std::vector<Eigen::Vector3d> points;



        // yanhao: use laser intensity
        std::vector<float> Intensity_oneScan;
        getPointsIntensity(*scan,Intensity_oneScan);

        float angle_increment = getAngleIncrement(*scan);   // laser angle increment in [rad]
        angle_increment = angle_increment * (180.0f/ M_PI);
        float camera_fov = 120.0f;   // assume the camera field of view is 120 degree

        // Yanhao debug
//        std::cout<< "Yanhao: Check: angle_increment --- \n" <<std::endl;
//        std::cout<< angle_increment <<std::endl;
        points = AutoGetLinePts(Points,angle_increment,false);     // this is only for one scan

        std::vector<Eigen::Vector3d> points_yh;
        std::vector<float> points_intensity_yh;
        AutoGetLinePts_useReflector(points_yh, points_intensity_yh, Points,Intensity_oneScan, camera_fov, angle_increment, true); // use intensity to get points from calibraiton board




        // yanhao intensity of points


        // yanhao
//        Eigen::Vector3d apoint;
//        for(int i=0; i<10;i++){
//            apoint = points.at(i);
//            std::cout<< "Yanhao: Check: point \n" <<std::endl;
//            std::cout<<"apoint["<<i<<"]="<< apoint <<std::endl;
//        }

      // 检测到了直线
      if(points.size() > 0)
      {
        // 在 camera 里找时间戳最近的一个 pose
        double min_dt = 10000;
        CamPose colsetTagPose;
        for (int i = 0; i < tagpose.size(); ++i) {
          CamPose tmp = tagpose.at(i);
          double t = fabs(tmp.timestamp - timestamp);
          if(t < min_dt)
          {
            min_dt = t;
            colsetTagPose = tmp;
          }
        }


        if(min_dt < 0.02)  // 20ms
        {

//          std::cout << "scan and tag time: "<<std::fixed<<std::setprecision(18)
//                    <<timestamp<<" "<<colsetTagPose.timestamp<<std::endl;
          /////////////////////////////////////////////////

          Eigen::Vector2d line;
          LineFittingCeres(points,line);
          std::vector<Eigen::Vector3d> points_on_line;

          // 激光所在直线不能垂直于某个轴
          double x_start(points.begin()->x()), x_end(points.end()->x());
          double y_start(points.begin()->y()), y_end(points.end()->y());
          if( fabs(x_end - x_start) > fabs(y_end - y_start) )
          {
            y_start = - (x_start * line(0) + 1) / line(1);
            y_end = - (x_end * line(0) + 1) / line(1);

          } else // 可能垂直于 x 轴，采用y值来计算 x
          {
            x_start = - (y_start * line(1) + 1) / line(0);
            x_end = - (y_end * line(1) + 1) / line(0);
          }

          points_on_line.push_back(Eigen::Vector3d(x_start,y_start,0));
          points_on_line.push_back(Eigen::Vector3d(x_end,y_end,0));

          Oberserve ob;
          ob.tagPose_Qca = colsetTagPose.qwc.inverse();          // tagPose_Qca　存的是世界在相机坐标系的坐标 camera (origin), 世界指版坐标系
          ob.tagPose_tca = -ob.tagPose_Qca.toRotationMatrix() * colsetTagPose.twc;
          ob.points = points;
          ob.points_on_line = points_on_line;


          // yanhao: remove unreasonable line fitting
          bool flag_push_this_ob = false;
          double lineFit_threshold = 10.0;

          Eigen::Vector3d points_on_line_1 = points_on_line[0];
          Eigen::Vector3d points_on_line_2 = points_on_line[1];

          if(std::fabs(points_on_line_1(0))<lineFit_threshold &  std::fabs(points_on_line_1(1))<lineFit_threshold &
                  std::fabs(points_on_line_2(0))<lineFit_threshold &  std::fabs(points_on_line_2(1))<lineFit_threshold ) { // check whether it is reasonable
            flag_push_this_ob = true;
          }


//          for(size_t k=0; k<points_on_line.size(); k++){
//            Eigen::Vector3d tmp_points_on_line = points_on_line[k];
//            if(tmp_points_on_line(0)<lineFit_threshold & tmp_points_on_line(1)<lineFit_threshold){ // check whether it is reasonable
//
//
//
//              std::cout<< "Yanhao: Check:----- \n" <<std::endl;
//              std::cout<< "before flag_push_this_ob = "<<flag_push_this_ob <<std::endl;
//              std::cout<< "tmp_points_on_line(0) = "<<tmp_points_on_line(0)<<"   " << "tmp_points_on_line(1) = "<<tmp_points_on_line(1) <<std::endl;
//              std::cout<< "after flag_push_this_ob = "<<flag_push_this_ob <<std::endl;
//
//              flag_push_this_ob = true;
//            }
//          }

          if(flag_push_this_ob){
            obs.push_back(ob);
//            flag_push_this_ob = false;
          }


        }

      }

    }
  }

    bool debug_yh = true;
    std::string saveobs_path = "/home/yanhao/LaserCameraCal_ws/data/";
    if(debug_yh){
        saveobs(obs,saveobs_path);
    }



  if(obs.size() < 5)
  {
    std::cout << "Valid Calibra Data Less"<<std::endl;
    bag_input.close();
    return 0;
  }
  std::cout <<"obs size: "<< obs.size() <<std::endl;

//    Eigen::Matrix4d Tcl0 = Eigen::Matrix4d::Identity();
    // Yanhao Chang: identity is not a good initial value
    Eigen::Matrix4d Tlc0 = Eigen::Matrix4d::Identity();
//    Tlc0(0,3) = 0.0; // -27.2/1000.0;
//    Tlc0(1,3) = 0.0; // -15.1/1000.0;
//    Tlc0(2,3) = 0.0; // 57.4/1000.0;
    Tlc0(0,3) =  -27.2/1000.0;   // x
    Tlc0(1,3) =  -15.1/1000.0;   // y
    Tlc0(2,3) =  57.4/1000.0;    // z

    Eigen::Matrix3d m;
    m << 0, 0, 1,
         -1, 0, 0,
         0, -1, 0;
    Tlc0.block<3,3>(0,0) = m;
    Eigen::Matrix4d Tcl = Tlc0.inverse();

    std::cout<< "Yanhao: Check:----- \n" <<std::endl;
    std::cout<< Tlc0 <<std::endl;

  //CamLaserCalibration(obs,Tcl, false);
  CamLaserCalibration(obs,Tcl, false, true);     // 建议用false true 或者 true true

//    std::cout << "\n----- Transform from laser to camera  Tcl is: -----\n"<<std::endl;
//    std::cout<< Tcl <<std::endl;

  std::cout << "\n----- Transform from Camera to Laser (origin) Tlc is: -----\n"<<std::endl;
  Eigen::Matrix4d Tlc = Tcl.inverse();
  std::cout<< Tlc <<std::endl;

  std::cout << "\n----- Transform from Camera to Laser, euler angles and translations are: -----\n"<<std::endl;
  Eigen::Matrix3d Rlc(Tlc.block(0,0,3,3));
  Eigen::Vector3d tlc(Tlc.block(0,3,3,1));
  EulerAngles rpy =  ToEulerAngles(Eigen::Quaterniond(Rlc));
  std::cout << "   roll(rad): "<<rpy.roll <<" pitch(rad): "<<rpy.pitch << " yaw(rad): "<<rpy.yaw<<"\n"
            << "or roll(deg): "<<rpy.roll * 180./M_PI <<" pitch(deg): "<<rpy.pitch* 180./M_PI  << " yaw(deg): "<<rpy.yaw * 180./M_PI <<"\n"
            << "       tx(m): "<<tlc.x() << "  ty(m): "<<tlc.y() << "   tz(m): "<<tlc.z()<<std::endl;

  // save to yaml file
  cv::Mat cvTlc;
  cv::eigen2cv(Tlc,cvTlc);
  std::string fn = savePath + "result.yaml";
  cv::FileStorage fs(fn, cv::FileStorage::WRITE);
  fs << "extrinsicTlc"<<cvTlc;
  cv::Mat cvrpy;
  cv::eigen2cv(Eigen::Vector3d(rpy.roll,rpy.pitch,rpy.yaw),cvrpy);
  cv::Mat cvtlc;
  cv::eigen2cv(tlc,cvtlc);
  fs << "RollPitchYaw"<<cvrpy;
  fs << "txtytz"<<cvtlc;
  fs.release();


  std::cout << "\n Result file : "<<fn<<std::endl;
  std::cout << "\n-------------- Calibration Code End --------------\n"<<std::endl;


  ros::spin();
}
