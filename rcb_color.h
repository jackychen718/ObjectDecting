#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "opencv2/opencv.hpp"
#include <thread>

using namespace cv;
using namespace std;

#define NUM_SAMPLES 200 //ÿ�����ص���������� 
#define SUBSAMPLE_FACTOR 8 //�Ӳ������� 
#define CHANNEL_THRO 8 // ����ȷ����Ԫ��ͨ������ֵ

#define TRUE 1
#define FALSE 0
#define u_char unsigned char
#define CHANNELS 3
#define MAXTAG 3    
#define THREAD_NUM (8)


typedef struct ce {
    u_char learnHigh[CHANNELS]; // ����Ԫ��ͨ���ķ�ֵ����(ѧϰ����) 
    u_char learnLow[CHANNELS]; // ����Ԫ��ͨ���ķ�ֵ���� 
    int t_last_update; // ����Ԫ���һ�θ��µ�ʱ��,ÿһ֡Ϊһ����λʱ��,���ڼ���stale 
    int stale; // ����Ԫ�������ʱ��,����ɾ���涨ʱ�䲻���µ���Ԫ,�����뱾 
} code_element; // ��Ԫ�����ݽṹ 

typedef struct code_book {
    code_element **cb;  // ��Ԫ�Ķ�άָ��,���Ϊָ����Ԫָ�������ָ��,ʹ�������Ԫʱ����Ҫ���ظ�����Ԫ,ֻ��Ҫ�򵥵�ָ�븳ֵ���� 
    int numEntries; // ���뱾����Ԫ����Ŀ 
    int t; // ���뱾���ڵ�ʱ��,һ֡Ϊһ��ʱ�䵥λ 
    int unmatch;    // û��ƥ��Ĵ���
} codeBook; // �뱾�����ݽṹ 

class Random_CB
{
public:
    Random_CB(void);
    ~Random_CB(void);

    void init(const Mat _image, int core); //��ʼ�� 
    void processFirstFrame(u_char *p, int width, int height, int numChannels, int core);
    void testAndUpdate(u_char *p, int width, int height, int numChannels, int core); //���� 
    Mat getMask(int core) { return m_mask[core]; };
    void deleteSamples() { delete [] cB; };
    int cvclearStaleEntries(codeBook &c);
    void postprocess(Mat frame, Mat mask, int width, int hight);

    std::thread *tasks[THREAD_NUM];
    
    cv::Mat frame;
    void executeThread(cv::Mat frame);
    void executeProc(); //�㷨ִ�������
    void MultiThread(const cv::Mat frame);
    void ThreadProc(void *arg);
    Mat maskAll;

private:
    codeBook *cB[THREAD_NUM];
    Mat m_mask[THREAD_NUM];
    bool flag[THREAD_NUM]; //��һ֡��ʶ
};

typedef struct _process_info
{
    int Index;
    int High_H;
    int Wigh_W;
    uchar *pSrcData;
} ProcessInfo;

