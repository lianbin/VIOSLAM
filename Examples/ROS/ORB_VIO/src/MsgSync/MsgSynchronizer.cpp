#include "MsgSynchronizer.h"
#include "../../../src/IMU/configparam.h"

namespace ORBVIO
{

MsgSynchronizer::MsgSynchronizer(const double& imagedelay):
    _imageMsgDelaySec(imagedelay), _status(NOTINIT),
    _dataUnsyncCnt(0)
{
    printf("image delay set as %.1fms\n",_imageMsgDelaySec*1000);//秒->毫秒
}

MsgSynchronizer::~MsgSynchronizer()
{

}


bool MsgSynchronizer::getRecentMsgs(sensor_msgs::ImageConstPtr &imgmsg, std::vector<sensor_msgs::ImuConstPtr> &vimumsgs)
{
    //unique_lock<mutex> lock2(_mutexIMUQueue);
    unique_lock<mutex> lock1(_mutexImageQueue);

    if(_status == NOTINIT || _status == INIT)
    {
        //ROS_INFO("synchronizer not inited");
        return false;
    }
    if(_imageMsgQueue.empty())
    {
        //ROS_INFO("no image stored in queue currently.");
        return false;
    }
    if(_imuMsgQueue.empty())
    {
        //ROS_WARN("no imu message stored, shouldn't");
        return false;
    }

    {
        sensor_msgs::ImageConstPtr imsg;
        sensor_msgs::ImuConstPtr bmsg;

        //
        imsg = _imageMsgQueue.back(); //最新的image
        bmsg = _imuMsgQueue.front();  //最老的imu

        // Check dis-continuity, tolerance 3 seconds
        //IMU数据超前image数据3秒以上，清空缓存数据
        if(imsg->header.stamp.toSec()-_imageMsgDelaySec + 3.0 < bmsg->header.stamp.toSec() )
        {
            ROS_ERROR("Data dis-continuity, > 3 seconds. Buffer cleared");
            clearMsgs();
            return false;
        }

        //
        imsg = _imageMsgQueue.front();//最老的image
        bmsg = _imuMsgQueue.back();   //最新的imu

        // Wait imu messages in case communication block
        if(imsg->header.stamp.toSec()-_imageMsgDelaySec > bmsg->header.stamp.toSec())
        {
            return false;
        }

        // Check dis-continuity, tolerance 3 seconds
        //这个判断如果满足，则一定满足上一个判断，所以找个判断应该没有用
        if(imsg->header.stamp.toSec()-_imageMsgDelaySec > bmsg->header.stamp.toSec() + 3.0)
        {
            ROS_ERROR("Data dis-continuity, > 3 seconds. Buffer cleared");
            clearMsgs();
            return false;
        }

        // Wait until the imu packages totolly com
        // 图像帧是一定小于10个的，因为buffer就开了2个
        if(_imageMsgQueue.size()<10 && _imuMsgQueue.size()<15
           && imsg->header.stamp.toSec()-_imageMsgDelaySec>bmsg->header.stamp.toSec() )
        {
            //ROS_WARN_STREAM("here return, last imu time "<<);
            return false;

        }

    }

    // get image message
    imgmsg = _imageMsgQueue.front();//取最老的图像
    _imageMsgQueue.pop();           

    // clear imu message vector, and push all imu messages whose timestamp is earlier than image message
    vimumsgs.clear();
    while(true)
    {
        // if no more imu messages, stop loop
        if(_imuMsgQueue.empty())
            break;

        // consider delay between image and imu serial
        sensor_msgs::ImuConstPtr tmpimumsg = _imuMsgQueue.front();
		//将图像帧之前的IMU数据push到vimumsgs中
        if(tmpimumsg->header.stamp.toSec() < imgmsg->header.stamp.toSec() - _imageMsgDelaySec)
        {
            // add to imu message vector
            vimumsgs.push_back(tmpimumsg);
            {
                unique_lock<mutex> lock(_mutexIMUQueue);
                _imuMsgQueue.pop();
            }

            _dataUnsyncCnt = 0;
        }
        else
        {
            if(_dataUnsyncCnt++>10)
            {
                _dataUnsyncCnt = 0;
                //_imuMsgQueue = std::queue<sensor_msgs::ImuConstPtr>();
                clearMsgs();
                ROS_ERROR("data unsynced many times, reset sync");
                return false;
            }
            // stop loop
            break;
        }
    }

    // the camera fps 20Hz, imu message 100Hz. so there should be not more than 5 imu messages between images
    if(vimumsgs.size()>10)//vimumsgs中存储的是两帧之间的IMU数据
        ROS_WARN("%lu imu messages between images, note",vimumsgs.size());
    if(vimumsgs.size()==0)
        ROS_ERROR("no imu message between images!");

    return true;
}

void MsgSynchronizer::addImuMsg(const sensor_msgs::ImuConstPtr &imumsg)
{
    unique_lock<mutex> lock(_mutexIMUQueue);

    if(_imageMsgDelaySec>=0) {
        _imuMsgQueue.push(imumsg);
        if(_status == NOTINIT)
        {
            _imuMsgTimeStart = imumsg->header.stamp;
            _status = INIT;
        }
    }
    else {
        // if there's no image messages, don't add image
        if(_status == NOTINIT)
            return;
        else if(_status == INIT)
        {
            // ignore all image messages with no imu messages between them
            // only add below images
            if(imumsg->header.stamp.toSec() + _imageMsgDelaySec > _imuMsgTimeStart.toSec())
            {
                _imuMsgQueue.push(imumsg);
                _status = NORMAL;
            }
        }
        else
        {
            // push message into queue
            _imuMsgQueue.push(imumsg);
        }
    }


}

//图像回调函数
void MsgSynchronizer::addImageMsg(const sensor_msgs::ImageConstPtr &imgmsg)
{
    unique_lock<mutex> lock(_mutexImageQueue);

    if(_imageMsgDelaySec >= 0) {//图像数据延后与IMU数据
        // if there's no imu messages, don't add image
        if(_status == NOTINIT)//图像数据延后的情况下，必然是先产生IMU数据，所以没有IMU数据，则图像数据也舍弃。
            return;
        else if(_status == INIT)
        {
            // ignore all image messages with no imu messages between them
            // only add below images
            //在获取到imu数据的时间点之前的图像，全都舍弃掉
            if(imgmsg->header.stamp.toSec() - _imageMsgDelaySec > _imuMsgTimeStart.toSec())
            {
                _imageMsgQueue.push(imgmsg);
                _status = NORMAL;
            }
        }
        else
        {
            // push message into queue
            _imageMsgQueue.push(imgmsg);
        }
    }
    else {  // start by image message ,delay <0 则图像数据先产生，应该根据图像数据进行初始化
        if(_status == NOTINIT)
        {
            _imuMsgTimeStart = imgmsg->header.stamp;
            _status = INIT;
        }
        else
        {   // no image data if there's no imu message
            _imageMsgQueue.push(imgmsg);
        }

    }

    if(ORB_SLAM2::ConfigParam::GetRealTimeFlag())
    {
        // Ignore earlier frames
        if(_imageMsgQueue.size()>2)//队列中只保留两个图像帧
            _imageMsgQueue.pop();
    }
}


void MsgSynchronizer::imageCallback(const sensor_msgs::ImageConstPtr& msg)
{
    addImageMsg(msg);
}

void MsgSynchronizer::imuCallback(const sensor_msgs::ImuConstPtr &msg)
{
    addImuMsg(msg);
}

void MsgSynchronizer::clearMsgs(void)
{
    _imuMsgQueue = std::queue<sensor_msgs::ImuConstPtr>();
    _imageMsgQueue = std::queue<sensor_msgs::ImageConstPtr>();
//    while(!_imageMsgQueue.empty())
//    {
//        _imageMsgQueue.pop();
//    }
//    while(!_imuMsgQueue.empty())
//    {
//        _imuMsgQueue.pop();
//    }
}

}
