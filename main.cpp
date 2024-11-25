#include "rcb_color.h"
#include <cstdio>

using namespace cv;
using namespace std;

int main(int argc, char* argv[])
{
    cv::Mat frame;
    VideoCapture capture;
    capture.open("video\\0005.avi");
	//capture.open("0005-120.avi");

    if (!capture.isOpened()){
        cout << "No camera or video input!\n" << endl;
        system("pause");
        return -1;
    }

    Random_CB RandomCB;
    while (capture.read(frame)) {
        if (frame.empty())
            break;

        RandomCB.frame = frame;
        RandomCB.executeProc();

        imshow("frame", frame);
       
        if (27 == cvWaitKey(1))
            break;
    }
    return 0;
}