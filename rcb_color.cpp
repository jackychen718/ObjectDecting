#include "rcb_color.h"

int c_xoff[9] = { -1, 0, 1, -1, 1, -1, 0, 1, 0 }; //x���ھӵ�
int c_yoff[9] = { -1, 0, 1, -1, 1, -1, 0, 1, 0 }; //y���ھӵ�

Random_CB::Random_CB(void)
{
    for(int i=0; i<THREAD_NUM; i++){
        flag[i] = true;
    } 
}

Random_CB::~Random_CB(void)
{
    for (size_t i = 0; i < THREAD_NUM; ++i) { 
    if (tasks[i]->joinable()) 
        tasks[i]->join(); 
        delete tasks[i]; 
    }
}

void Random_CB::init(const Mat _image, int core)
{
    cB[core] = (codeBook*)calloc(_image.rows*_image.cols, sizeof(codeBook));
    for (int i = 0; i < _image.rows*_image.cols; i++){
        cB[core][i].numEntries = 0;
        cB[core][i].unmatch = 0;
    }
    m_mask[core] = Mat::zeros(_image.size(), CV_8UC1);
}

void Random_CB::processFirstFrame(u_char *p, int width, int height, int numChannels, int core)
{
    RNG rng;
    int row, col;
    unsigned cbBounds[CHANNELS];
    int high[3], low[3];
    int matchChannel;

    for (int i = 0; i < numChannels; i++){
        cbBounds[i] = CHANNEL_THRO; // ����ȷ����Ԫ��ͨ���ķ�ֵ 
    }
    FILE *found;
    found = fopen("codeelement_num.txt", "w");

    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            for (int k = 0; k < NUM_SAMPLES; k++){
                int ii;
                int random = rng.uniform(0, 9);

                row = i + c_yoff[random];
                if (row < 0)
                    row = 0;
                if (row >= height)
                    row = height - 1;

                col = j + c_xoff[random];
                if (col < 0)
                    col = 0;
                if (col >= width)
                    col = width - 1;

                if (cB[core][i*width + j].numEntries == 0){
                    cB[core][i*width + j].t = 0;
                }

                cB[core][i*width + j].t += 1;
                for (int n = 0; n < numChannels; n++) {
                    high[n] = *(p + CHANNELS*row*width + CHANNELS * col + n) + cbBounds[n];
                    if (high[n] > 255) high[n] = 255;
                    low[n] = *(p + CHANNELS * row*width + CHANNELS * col + n) - cbBounds[n];
                    if (low[n] < 0) low[n] = 0;
                }
                for (ii = 0; ii < cB[core][i*width + j].numEntries; ii++) {
                    matchChannel = 0;
                    for (int n = 0; n < numChannels; n++) {
                        if ((cB[core][i*width + j].cb[ii]->learnLow[n] <= *(p + CHANNELS * row*width + CHANNELS * col + n)) && (*(p + CHANNELS * row*width + CHANNELS * col + n) <= cB[core][i*width + j].cb[ii]->learnHigh[n])){ //Found an entry for this channel 
                            matchChannel++;
                        }
                    }
                    if (matchChannel == numChannels) {
                        cB[core][i*width + j].cb[ii]->t_last_update = cB[core][i*width + j].t;
                        for (int n = 0; n < numChannels; n++){
                            if (cB[core][i*width + j].cb[ii]->learnHigh[n] < high[n])
                                cB[core][i*width + j].cb[ii]->learnHigh[n] = high[n];
                            else if (cB[core][i*width + j].cb[ii]->learnLow[n] > low[n])
                                cB[core][i*width + j].cb[ii]->learnLow[n] = low[n];
                        }
                        break;
                    }
                }

                code_element **foo[THREAD_NUM];
                if (ii == cB[core][i*width + j].numEntries){
                    foo[core] = new code_element*[cB[core][i*width + j].numEntries + 1];
                    for (int jj = 0; jj < cB[core][i*width + j].numEntries; jj++)
                        foo[core][jj] = cB[core][i*width + j].cb[jj];

                    foo[core][cB[core][i*width + j].numEntries] = new code_element;
                    if (cB[core][i*width + j].numEntries) delete[] cB[core][i*width + j].cb;
                    cB[core][i*width + j].cb = foo[core];
                    for (int n = 0; n < numChannels; n++){
                        cB[core][i*width + j].cb[cB[core][i*width + j].numEntries]->learnHigh[n] = high[n];
                        cB[core][i*width + j].cb[cB[core][i*width + j].numEntries]->learnLow[n] = low[n];
                    }
                    cB[core][i*width + j].cb[cB[core][i*width + j].numEntries]->t_last_update = cB[core][i*width + j].t;
                    cB[core][i*width + j].cb[cB[core][i*width + j].numEntries]->stale = 0;
                    cB[core][i*width + j].numEntries += 1;
                }
                for (int s = 0; s < cB[core][i*width + j].numEntries; s++){
                    int negRun = cB[core][i*width + j].t - cB[core][i*width + j].cb[s]->t_last_update;
                    if (cB[core][i*width + j].cb[s]->stale < negRun)
                        cB[core][i*width + j].cb[s]->stale = negRun;
                }
            }
            cvclearStaleEntries(cB[core][i*width + j]);
            fprintf(found, "%d\n", cB[core][i*width + j].numEntries);
        }
    }
    fclose(found);
}


void Random_CB::testAndUpdate(u_char *p, int width, int height, int numChannels, int core)
{
    RNG rng;
    unsigned cbBounds[CHANNELS];
    int high[3], low[3];
    int matchChannel;

    for (int i = 0; i < numChannels; i++){
        cbBounds[i] = CHANNEL_THRO; // ����ȷ����Ԫ��ͨ���ķ�ֵ 
    }
    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            int ii;
            for (int n = 0; n<numChannels;n++){
                high[n] = *(p + CHANNELS * i*width + CHANNELS * j + n) + cbBounds[n];// *(p+n) �� p[n] ����ȼ�,������*(p+n) �ٶȸ��� 
                if (high[n] > 255) high[n] = 255;
                low[n] = *(p + CHANNELS * i*width + CHANNELS * j + n) - cbBounds[n];
                if (low[n] < 0) low[n] = 0;
                // ��p ��ָ����ͨ������,�Ӽ�cbBonds����ֵ,��Ϊ�����ط�ֵ�������� 
            }
            cB[core][i*width + j].t += 1;
            for (ii = 0; ii < cB[core][i*width + j].numEntries; ii++){
                matchChannel = 0;
                for (int n = 0; n < numChannels; n++){
                    //����ÿ��ͨ�� 
                    if ((cB[core][i*width + j].cb[ii]->learnLow[n] <= *(p + CHANNELS * i*width + CHANNELS * j + n)) && (*(p + CHANNELS * i*width + CHANNELS * j + n) <= cB[core][i*width + j].cb[ii]->learnHigh[n])){ //Found an entry for this channel                                                                                                                                                                                   // ���p ����ͨ�������ڸ���Ԫ��ֵ������֮�� 
                        matchChannel++;
                    }
                }
                if (matchChannel == numChannels){// ���p ����ͨ�������ڸ���Ԫ��ֵ������֮�� 
                    m_mask[core].at<uchar>(i, j) = 0;
                    int random = rng.uniform(0, SUBSAMPLE_FACTOR);
                    if (random == 0){
                        random = rng.uniform(0, cB[core][i*width + j].numEntries);
                        cB[core][i*width + j].cb[random]->t_last_update = cB[core][i*width + j].t;
                        // ���¸���Ԫʱ��Ϊ��ǰʱ��
                        for (int n = 0; n < numChannels; n++){
                            if (cB[core][i*width + j].cb[random]->learnHigh[n] < high[n])
                                cB[core][i*width + j].cb[random]->learnHigh[n] = high[n];
                            else if (cB[core][i*width + j].cb[random]->learnLow[n] > low[n])
                                cB[core][i*width + j].cb[random]->learnLow[n] = low[n];
                        }

                    }
                    random = rng.uniform(0, SUBSAMPLE_FACTOR);
                    if (0 == random){
                        int row, col;
                        int random = rng.uniform(0, 9);

                        row = i + c_yoff[random];
                        if (row < 0)row = 0;
                        if (row >= height)
                            row = height - 1;

                        col = j + c_xoff[random];
                        if (col < 0)
                            col = 0;
                        if (col >= width)
                            col = width - 1;

                        random = rng.uniform(0, cB[core][row*width + col].numEntries);
                        cB[core][row*width + col].cb[random]->t_last_update = cB[core][i*width + j].t;
                        // ���¸���Ԫʱ��Ϊ��ǰʱ��
                        for (int n = 0; n < numChannels; n++) {
                            if (cB[core][row*width + col].cb[random]->learnHigh[n] < high[n])
                                cB[core][row*width + col].cb[random]->learnHigh[n] = high[n];
                            if (cB[core][row * width + col].cb[random]->learnLow[n] > low[n])
                                cB[core][row*width + col].cb[random]->learnLow[n] = low[n];
                        }
                    }
                    break;
                }
            }
            
            if (ii == cB[core][i*width + j].numEntries){
                cB[core][i*width + j].unmatch += 1;
                m_mask[core].at<uchar>(i, j) = 255;
                if (cB[core][i*width + j].unmatch > 50){
                    int random = rng.uniform(0, SUBSAMPLE_FACTOR);// �޸ĺ�ģ�Ӧ����������ʸ��µ�ǰֵ��ģ��
                    if (0 == random){
                        random = rng.uniform(0, cB[core][i*width + j].numEntries);
                        cB[core][i*width + j].cb[random]->t_last_update = cB[core][i*width + j].t;

                        for (int n = 0; n < numChannels; n++){
                            if (cB[core][i*width + j].cb[random]->learnHigh[n] < high[n])
                                cB[core][i*width + j].cb[random]->learnHigh[n] = high[n];
                            else if (cB[core][i*width + j].cb[random]->learnLow[n] > low[n])
                                cB[core][i*width + j].cb[random]->learnLow[n] = low[n];
                        }
                        cB[core][i*width + j].unmatch = 0; //���º�����Ϊ0
                    }
                }
            }
            
            // ���¸���ʱ��
            for (int s = 0; s < cB[core][i*width + j].numEntries; s++){
                // This garbage is to track which codebook entries are going stale
                int negRun = cB[core][i*width + j].t - cB[core][i*width + j].cb[s]->t_last_update;
                // �������Ԫ�Ĳ�����ʱ�� 
                if (cB[core][i*width + j].cb[s]->stale < negRun)
                    cB[core][i*width + j].cb[s]->stale = negRun;
            }
            //cvclearStaleEntries(cB[core][i*width + j]);
        }
    }
}

int Random_CB::cvclearStaleEntries(codeBook &c)
{
    int staleThresh = floor(c.t / c.numEntries) - 2; // �趨ˢ��ʱ�� 
    if (staleThresh < 0) staleThresh = 1;

    int *keep = new int[c.numEntries]; // ����һ��������� 
    int keepCnt = 0; // ��¼��ɾ����Ԫ��Ŀ//SEE WHICH CODEBOOK ENTRIES ARE TOO STALE 
    for (int i = 0; i<c.numEntries;i++){ // �����뱾��ÿ����Ԫ 
        if (c.cb[i]->stale > staleThresh) {
            // ����Ԫ�еĲ�����ʱ������趨��ˢ��ʱ��,����Ϊɾ�� 
            keep[i] = 0; //Mark for destruction 
            delete c.cb[i];
        }
        else{
            keep[i] = 1; //Mark to keep 
            keepCnt += 1;
        }
    }// KEEP ON LY THE GOOD 
    c.t = 0; //Full reset on stale tracking 
    // �뱾ʱ������ 
    code_element **foo = new code_element*[keepCnt];
    // �����СΪkeepCnt ����Ԫָ������ 
    int k = 0;
    for (int ii = 0; ii < c.numEntries; ii++){
        if (keep[ii]){
            foo[k] = c.cb[ii];
            foo[k]->stale = 0; //We have to refresh these entries for next clearStale 
            foo[k]->t_last_update = 0;
            k++;
        }
    }//CLEAN UP delete[] keep;
    delete[] c.cb;
    delete[] keep;
    c.cb = foo;
    // ��foo ͷָ���ַ����c.cb 
    int numCleared = c.numEntries - keepCnt;
    // ���������Ԫ���� 
    c.numEntries = keepCnt;
    // ʣ�����Ԫ��ַ 
    return numCleared;
}


void Random_CB::postprocess(Mat frame, Mat mask, int width, int hight)
{
    Mat mask_dilate;
    Mat kernel = getStructuringElement(MORPH_RECT, Size(19, 19), Point(-1, -1));
    dilate(mask, mask_dilate, kernel, Point(-1, -1), 1);
    
    vector<vector<Point>> contours;
    vector<Vec4i> hireachy;
    findContours(mask_dilate, contours, hireachy, RETR_TREE, CHAIN_APPROX_SIMPLE, Point());

    Mat overlay;
    frame.copyTo(overlay);

    for (size_t t = 0; t < contours.size() && t < MAXTAG; t++) {
        drawContours(overlay, contours, t, Scalar(0, 97, 255), -1, 8, Mat(), 0, Point());
    }
    cv::addWeighted(overlay, 0.4, frame, 0.6, 0, frame);
}

void Random_CB::ThreadProc(void *arg)
{
    ProcessInfo *Info = (ProcessInfo *)arg;
    cv::Mat frame = Mat(Info->High_H/THREAD_NUM, Info->Wigh_W, CV_8UC1, Info->pSrcData);

    //cout << " [ThreadProc] Wide=" <<frame.cols<<" High="<<frame.rows<< endl;

    if (flag[Info->Index]) {
        init(frame, Info->Index);
        processFirstFrame(frame.data, frame.cols, frame.rows, CHANNELS, Info->Index);
        flag[Info->Index] = false;
    } else {
        testAndUpdate(frame.data, frame.cols, frame.rows, CHANNELS, Info->Index);
    }
    return;
}

void Random_CB::MultiThread(const cv::Mat frame)
{
     //cout << " [MultiThread] Wide=" <<frame.cols<<" High="<<frame.rows<<"channels="<<frame.channels()<<endl;
    if (false == frame.isContinuous()){
        cout << " [error] frame is not Continuous" << endl;
        return;
    }

    int block_size = frame.cols*frame.rows*frame.channels() / THREAD_NUM;
    int Offset = 0;
    ProcessInfo InfoArray[THREAD_NUM];
    
    for (int i = 0; i < THREAD_NUM; i++) {
        Offset = block_size * i;
        InfoArray[i].Index = i;
        InfoArray[i].Wigh_W = frame.cols;
        InfoArray[i].High_H = frame.rows;
        InfoArray[i].pSrcData = frame.data + Offset;
        tasks[i] = new std::thread(&Random_CB::ThreadProc, this, &InfoArray[i]);
    }

    for (int i = 0; i < THREAD_NUM; i++) {
        tasks[i]->join();
    }
    return;
}

void Random_CB::executeProc()
{
    cv::Mat garyFrame;

    if (1 == frame.channels()) {
        cvtColor(frame, frame, CV_GRAY2RGB);
    }
    
    if (1 == CHANNELS) {
        cvtColor(frame, garyFrame, CV_RGB2GRAY);
    } else if (3 == CHANNELS){
        cvtColor(frame, garyFrame, CV_BGR2YCrCb);
    }

    MultiThread(garyFrame);

    int Offset = 0;
    int block_size = frame.cols*frame.rows / THREAD_NUM;
    cv::Mat mask = Mat::zeros(frame.size(), CV_8UC1);
    
    for (int i=0; i<THREAD_NUM; i++){
        Offset = block_size * i;
        memcpy(mask.data + Offset, getMask(i).data, block_size);
    }

    morphologyEx(mask, mask, MORPH_OPEN, Mat());
    maskAll = mask;
    postprocess(frame, mask, frame.cols, frame.rows);

}

void Random_CB::executeThread(cv::Mat frame)
{
    this->frame = frame; 
    std::thread t(&Random_CB::executeProc, this);
    t.join();
}