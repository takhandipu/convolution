/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   fusion.cpp
 * Author: takh
 *
 * Created on February 5, 2018, 8:40 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

#define DEBUG 0

using namespace std;

class Fusion
{
private:
    int W1, H1, D1;
    int K, F, S, P;
    int W2, H2, D2;
    int W3, H3, D3;
    int ***input;
    int ****filter1;
    int *bias1;
    int ***output1;
    bool ***flag;
    int ****filter2;
    int *bias2;
    int ***output2;
public:
    Fusion(int W1, int H1, int D1, int K, int F, int S, int P)
    {
        this->W1 = W1;
        this->H1 = H1;
        this->D1 = D1;
        this->K = K;
        this->F = F;
        this->S = S;
        this->P = P;
        this->W2 = ((W1 - F + (2*P))/S) + 1;
        this->H2 = ((H1 - F + (2*P))/S) + 1;
        this->D2 = K;
        this->W3 = ((W2 - F + (2*P))/S) + 1;
        this->H3 = ((H2 - F + (2*P))/S) + 1;
        this->D3 = K;
        if(DEBUG)
        {
            cout<<"W1: "<<W1<<endl;
            cout<<"H1: "<<H1<<endl;
            cout<<"D1: "<<D1<<endl;
            cout<<"K: "<<K<<endl;
            cout<<"F: "<<F<<endl;
            cout<<"S: "<<S<<endl;
            cout<<"P: "<<P<<endl;
            cout<<"W2: "<<W2<<endl;
            cout<<"H2: "<<H2<<endl;
            cout<<"D2: "<<D2<<endl;
            cout<<"W3: "<<W3<<endl;
            cout<<"H3: "<<H3<<endl;
            cout<<"D3: "<<D3<<endl;
        }
        input = new int**[D1];
        for(int i = 0; i < D1; i++)
        {
            input[i] = new int*[H1+2*P];
            for(int j = 0; j < H1+2*P; j++)
            {
                input[i][j] = new int[W1+2*P];
                for(int k = 0; k <W1+2*P; k++)
                {
                    if(P<=k && k<W1-1+2*P && P<=j && j<H1-1+2*P){
                        input[i][j][k] = rand() % 10;
                    }
                    else {
                        input[i][j][k] = 0;
                    }
                }
            }
        }
        if(DEBUG)
        {
            for(int i = 0; i < D1; i++)
            {
                cout<<"input[][]["<<i<<"]"<<endl;
                for(int j = 0; j < H1+2*P; j++)
                {
                    for(int k = 0; k <W1+2*P; k++)
                    {
                        cout<<input[i][j][k]<<" ";
                    }
                    cout<<endl;
                }
            }
        }
        output1 = new int **[D2];
        for(int i = 0; i < D2; i++)
        {
            output1[i] = new int *[H2+2*P];
            for(int j = 0; j < H2 + 2*P; j++)
            {
                output1[i][j] = new int[W2+2*P];
                for(int k = 0; k < W2+2*P; k++)
                {
                    output1[i][j][k] = 0;
                }
            }
        }
        flag = new bool **[D2];
        for(int i = 0; i < D2; i++)
        {
            flag[i] = new bool *[H2+2*P];
            for(int j = 0; j < H2 + 2*P; j++)
            {
                flag[i][j] = new bool[W2+2*P];
                for(int k = 0; k < W2+2*P; k++)
                {
                    flag[i][j][k] = false;
                    if(P<=k && k<W2+P && P<=j && j<H2+P){
                        flag[i][j][k] = false;
                    }
                    else {
                        flag[i][j][k] = true;
                    }
                }
            }
        }
        filter1 = new int ***[K];
        for(int i = 0; i < K; i++)
        {
            filter1[i] = new int **[D1];
            for(int j = 0; j < D1; j++)
            {
                filter1[i][j] = new int *[F];
                for(int k = 0; k < F; k++)
                {
                    filter1[i][j][k] = new int[F];
                    for(int l = 0; l < F; l++)
                    {
                        filter1[i][j][k][l] = rand() % 5;
                    }
                }
            }
        }
        if(DEBUG)
        {
            for(int i = 0; i < K; i++)
            {
                for(int j = 0; j < D1; j++)
                {
                    cout<<"w"<<i<<"[][]["<<j<<"]\n";
                    for(int k = 0; k < F; k++)
                    {
                        for(int l = 0; l < F; l++)
                        {
                            cout<<filter1[i][j][k][l]<<" ";
                        }
                        cout<<endl;
                    }
                }
            }
        }
        bias1 = new int[K];
        for(int i = 0; i < K; i++)
        {
            bias1[i] = rand()%3;
        }
        if(DEBUG)
        {
            cout<<"Bias"<<endl;
            for(int i = 0; i < K; i++)
            {
                cout << bias1[i]<<" ";
            }
            cout<<endl;
        }
        output2 = new int **[D3];
        for(int i = 0; i < D3; i++)
        {
            output2[i] = new int *[H3+2*P];
            for(int j = 0; j < H3 + 2*P; j++)
            {
                output2[i][j] = new int[W3+2*P];
                for(int k = 0; k < W3+2*P; k++)
                {
                    output2[i][j][k] = 0;
                }
            }
        }
        filter2 = new int ***[K];
        for(int i = 0; i < K; i++)
        {
            filter2[i] = new int **[D2];
            for(int j = 0; j < D2; j++)
            {
                filter2[i][j] = new int *[F];
                for(int k = 0; k < F; k++)
                {
                    filter2[i][j][k] = new int[F];
                    for(int l = 0; l < F; l++)
                    {
                        filter2[i][j][k][l] = rand() % 5;
                    }
                }
            }
        }
        if(DEBUG)
        {
            for(int i = 0; i < K; i++)
            {
                for(int j = 0; j < D2; j++)
                {
                    cout<<"w"<<i<<"[][]["<<j<<"]\n";
                    for(int k = 0; k < F; k++)
                    {
                        for(int l = 0; l < F; l++)
                        {
                            cout<<filter2[i][j][k][l]<<" ";
                        }
                        cout<<endl;
                    }
                }
            }
        }
        bias2 = new int[K];
        for(int i = 0; i < K; i++)
        {
            bias2[i] = rand()%3;
        }
        if(DEBUG)
        {
            cout<<"Bias"<<endl;
            for(int i = 0; i < K; i++)
            {
                cout << bias2[i]<<" ";
            }
            cout<<endl;
        }
    }
    int sum1(int w2, int h2, int d2)
    {
        int retValue = bias1[d2];
        for(int i = 0; i < D1; i++)
        {
            for(int j = 0; j <F; j++)
            {
                for(int k = 0; k < F; k++)
                {
                    if(DEBUG)cout<<i<<j+S*h2<<k+S*w2<<i<<j<<k<<" ";
                    //cout<<"Input,"<<i<<","<<j+S*h2<<","<<k+S*w2<<endl;
                    retValue += input[i][j+S*h2][k+S*w2] * filter1[d2][i][j][k];
                }
            }
        }
        if(DEBUG)cout<<endl;
        return retValue;
    }
    /*void run1()
    {
        //start convolution
        for(int i = 0; i < D2; i++)
        {
            if(DEBUG)cout<<"output1[][]["<<i<<"]\n";
            for(int j = 0; j < H2; j++)
            {
                for(int k = 0; k < W2; k++)
                {
                    output1[i][j+P][k+P] = sum1(k,j,i);
                    flag[i][j+P][k+P]=true;
                    if(DEBUG)cout<<output1[i][j+P][k+P]<<" ";
                }
                if(DEBUG)cout<<"\n";
            }
        }
    }*/
    int sum2(int w2, int h2, int d2)
    {
        int retValue = bias2[d2];
        for(int i = 0; i < D2; i++)
        {
            for(int j = 0; j <F; j++)
            {
                for(int k = 0; k < F; k++)
                {
                    /*if(!flag[i][j+S*h2][k+S*w2])
                    {
                        output1[i][j+S*h2][k+S*w2] = sum1(k+S*w2-P,j+S*h2-P,i);
                        //cout<<"1,"<<i<<","<<j+S*h2<<","<<k+S*w2<<endl;
                        flag[i][j+S*h2][k+S*w2]=true;
                    }*/
                    //if(DEBUG)cout<<i<<j+S*h2<<k+S*w2<<i<<j<<k<<" ";
                    //retValue += output1[i][j+S*h2][k+S*w2] * filter2[d2][i][j][k];
                    if(!flag[i][j+S*h2][k+S*w2])retValue += sum1(k+S*w2-P,j+S*h2-P,i) * filter2[d2][i][j][k];
                }
            }
        }
        if(DEBUG)cout<<endl;
        return retValue;
    }
    void run()
    {
        //start convolution
        /*for(int i = 0; i < D3; i++)
        {
            if(DEBUG)cout<<"output2[][]"<<i<<"]\n";
            for(int j = 0; j < H3; j++)
            {
                for(int k = 0; k < W3; k++)
                {
                    output2[i][j+P][k+P] = sum2(k,j,i);
                    //cout<<"2,"<<i<<","<<j+P<<","<<k+P<<endl;
                    if(DEBUG)cout<<output2[i][j+P][k+P]<<" ";
                }
                if(DEBUG)cout<<"\n";
            }
        }*/
        for(int i = 0; i < D2; i++)
        {
            if(DEBUG)cout<<"output[][]["<<i<<"]\n";
            for(int j = 0; j < H2; j++)
            {
                for(int k = 0; k < W2; k++)
                {
                    //cout<<"Computing, "<<i<<","<<j+P<<","<<k+P<<endl;
                    //cout<<"Accessing,"<<endl;
                    int val = sum1(k,j,i), conH, conW;
                    //if(DEBUG)cout<<output1[i][j+P][k+P]<<" ";
                    //cout<<"Layer1,"<<i<<","<<j+P<<","<<k+P<<endl;
                    for(int ii = 0; ii < D3; ii++)
                    {
                        for(int jj = -F+1; jj <1; jj+=S)
                        {
                            for(int kk = -F+1; kk <1; kk+=S)
                            {
                                conH = jj+j+P;
                                conW = kk+k+P;
                                //cout<<"Layer2,"<<ii<<","<<conH<<","<<conW<<endl;
                                //cout<<"Filter2,"<<ii<<","<<i<<","<<-jj<<","<<-kk<<endl;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=val*filter2[ii][i][-jj][-kk];
                                    //cout<<"Layer2,"<<ii<<","<<conH<<","<<conW<<endl;
                                    //cout<<"Filter2,"<<ii<<","<<i<<","<<-jj<<","<<-kk<<endl;
                                }
                            }
                        }
                    }
                }
                if(DEBUG)cout<<"\n";
            }
        }
    }
};

/*
 * 
 */
int main(int argc, char** argv) {
    ifstream inFile("input.txt");
    int W1, H1, D1, K, F, S, P;
    inFile >> W1 >> H1 >> D1;//5 5 3
    inFile >> K >> F >> S >> P;//2 3 2 1
    Fusion fusion(W1, H1, D1, K, F, S, P);
    struct timeval tv;
    gettimeofday(&tv,NULL);
    long start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    fusion.run();
    gettimeofday(&tv,NULL);
    long end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    long elapsed = end - start;
    cout<<"s, "<<(1.0 * elapsed)/1000000<<endl;
    return 0;
    return 0;
}

