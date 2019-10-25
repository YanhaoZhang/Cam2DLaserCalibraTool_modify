std::vector< Eigen::Vector3d > AutoGetLinePts(const std::vector<Eigen::Vector3d> points, bool debug)
{
//    cv::Mat img(img_w, img_w, CV_8UC1, cv::Scalar::all(0));
    cv::Mat img(img_w, img_w, CV_8UC3, cv::Scalar(0,0,0));

    for (auto pt: points) {
        int col = (int)(pt.x() / z * focal + img_w/2);
        int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

        if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
            continue;

        cv::Vec3b color_value(255,0,0);
        img.at<cv::Vec3b>(row, col) = color_value;
//        img.at<uchar>(row, col) = 255;
    }

    /// detect and get line
    // 直接从每一帧激光的正前方开始搜索一定距离范围内的符合平面标定板形状的激光线段
    int n = points.size();
    int id = n/2;
//        std::cout << points.at(id).transpose() <<" "<<points.at(id+1).transpose() <<std::endl;
    // 假设每个激光点之间的夹角为0.3deg,
    // step 1: 如果有激光标定板，那么激光标定板必须出现在视野的正前方 120 deg 范围内(通常相机视野也只有 120 deg)，也就是左右各 60deg.
    int delta = 80/0.3;

    int id_left = std::min( id + delta, n-1);
    int id_right = std::max( id - delta , 0);

    int dist_left = points.at(id_left).norm();
    int dist_right = points.at(id_right).norm();

    // 逻辑别搞复杂了。
    std::vector<LineSeg> segs;
    double dist_thre = 0.05;
    int skip = 3;
    int currentPt = id_right;
    int nextPt = currentPt + skip;
    bool newSeg = true;
    LineSeg seg;
    for (int i = id_right; i < id_left - skip; i += skip) {

        if(newSeg)
        {
            seg.id_start = currentPt;
            seg.id_end = nextPt;
            newSeg = false;
        }

        double d1 = points.at(currentPt).head(2).norm();
        double d2 = points.at(nextPt).head(2).norm();
        double range_max = 100;
        if(d1 < range_max && d2 < range_max)    // 有效数据,  激光小于 100 m
        {
            if(fabs(d1-d2) < dist_thre)  //  8cm
            {
                seg.id_end = nextPt;

            } else
            {
                newSeg = true;
                Eigen::Vector3d dist = points.at(seg.id_start) - points.at(seg.id_end);
                if(dist.head(2).norm() > 0.2
                   && points.at(seg.id_start).head(2).norm() < 2
                   && points.at(seg.id_end).head(2).norm() < 2
                   && seg.id_end-seg.id_start > 50   )  // 至少长于 20 cm, 标定板不能距离激光超过2m, 标定板上的激光点肯定多余 50 个
                {
                    seg.dist = dist.head(2).norm();
                    segs.push_back(seg);
                }
            }

            currentPt = nextPt;        // 下一个点
            nextPt = nextPt + skip;
        } else
        {
            if(d1 > range_max)
            {
                currentPt = nextPt;        // 下一个点
            }
            nextPt = nextPt + skip;
        }
    }


    // 对 segs 的边界进行扩充
    for (int i = 0; i < segs.size(); ++i) {
        LineSeg tmp = segs.at(i);

        for (int j = 1; j < 4; ++j)
        {
            int boundaryPt = tmp.id_end + j;
            double d1 = points.at(tmp.id_end).head(2).norm();
            double d2 = points.at(boundaryPt).head(2).norm();
            if(fabs(d1-d2) < dist_thre)  //  8cm
            {
                segs.at(i).id_end = boundaryPt;
            }
        }

        for (int j = -1; j > -4; --j)
        {
            int boundaryPt = tmp.id_start + j;
            double d1 = points.at(tmp.id_start).head(2).norm();
            double d2 = points.at(boundaryPt).head(2).norm();
            if(fabs(d1-d2) < dist_thre)  //  8cm
            {
                segs.at(i).id_start = boundaryPt;
            }
        }

    }

    std::vector<Eigen::Vector3d> ptsLine;
    if(segs.size() > 0)
    {
        LineSeg bestLine;
        int maxpts = -1;
        for (int i = 0; i < segs.size(); ++i) {
            LineSeg tmp = segs.at(i);
            int cnt = tmp.id_end - tmp.id_start;
            if(cnt > maxpts)
            {
                bestLine = tmp;
                maxpts = cnt;
            }
        }

        for (int i = bestLine.id_start; i < bestLine.id_end+1; ++i) {
            ptsLine.push_back(points.at(i));
        }
    }

    if(debug)
    {
        for (int j = 0; j < segs.size(); ++j) {
            LineSeg tmp = segs.at(j);
            for (int i = tmp.id_start; i < tmp.id_end; ++i)
            {
                Eigen::Vector3d pt = points.at(i);
                int col = (int)(pt.x() / z * focal + img_w/2);
                int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

                if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
                    continue;

                cv::Vec3b color_value(0,255,0);
                img.at<cv::Vec3b>(row, col) = color_value;
            }
        }

        for (int j = 0; j < ptsLine.size(); ++j) {

            Eigen::Vector3d pt = ptsLine.at(j);
            int col = (int)(pt.x() / z * focal + img_w/2);
            int row = (int)(- pt.y() / z * focal + img_w/2);  // -Y/Z 加了一个负号, 是为了抵消针孔投影时的倒影效果

            if(col > img_w-1 || col< 0 || row > img_w-1 || row < 0)
                continue;

            cv::Vec3b color_value(0,0,200);
            img.at<cv::Vec3b>(row, col) = color_value;

        }

        cv::putText(img, "Detecting the Laser Points on the calibra planar!",
                    cvPoint(5,30), cv::FONT_HERSHEY_COMPLEX_SMALL, 0.7, cvScalar(255,255,255), 1, CV_AA);
        cv::imshow("ScanPoint",img);
        cv::waitKey(10);
    }

    return ptsLine;

}
