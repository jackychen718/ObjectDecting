#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "opencv2/opencv.hpp"
#include <thread>

using namespace cv;
using namespace std;

#define NUM_SAMPLES 200 //每个像素点的样本个数 
#define SUBSAMPLE_FACTOR 8 //子采样概率 
#define CHANNEL_THRO 8 // 用于确定码元各通道的阈值

#define TRUE 1
#define FALSE 0
#define u_char unsigned char
#define CHANNELS 3
#define MAXTAG 3    
#define THREAD_NUM (8)


typedef struct ce {
    u_char learnHigh[CHANNELS]; // 此码元各通道的阀值上限(学习界限) 
    u_char learnLow[CHANNELS]; // 此码元各通道的阀值下限 
    int t_last_update; // 此码元最后一次更新的时间,每一帧为一个单位时间,用于计算stale 
    int stale; // 此码元最长不更新时间,用于删除规定时间不更新的码元,精简码本 
} code_element; // 码元的数据结构 

typedef struct code_book {
    code_element **cb;  // 码元的二维指针,理解为指向码元指针数组的指针,使得添加码元时不需要来回复制码元,只需要简单的指针赋值即可 
    int numEntries; // 此码本中码元的数目 
    int t; // 此码本现在的时间,一帧为一个时间单位 
    int unmatch;    // 没有匹配的次数
} codeBook; // 码本的数据结构 

class Random_CB
{
public:
    Random_CB(void);
    ~Random_CB(void);

    void init(const Mat _image, int core); //初始化 
    void processFirstFrame(u_char *p, int width, int height, int numChannels, int core);
    void testAndUpdate(u_char *p, int width, int height, int numChannels, int core); //更新 
    Mat getMask(int core) { return m_mask[core]; };
    void deleteSamples() { delete [] cB; };
    int cvclearStaleEntries(codeBook &c);
    void postprocess(Mat frame, Mat mask, int width, int hight);

    std::thread *tasks[THREAD_NUM];
    
    cv::Mat frame;
    void executeThread(cv::Mat frame);
    void executeProc(); //算法执行总入口
    void MultiThread(const cv::Mat frame);
    void ThreadProc(void *arg);
    Mat maskAll;

private:
    codeBook *cB[THREAD_NUM];
    Mat m_mask[THREAD_NUM];
    bool flag[THREAD_NUM]; //第一帧标识
};

typedef struct _process_info
{
    int Index;
    int High_H;
    int Wigh_W;
    uchar *pSrcData;
} ProcessInfo;

