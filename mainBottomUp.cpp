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
#include <algorithm>

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
                        input[i][j][k] = (i+j+k) % 10;
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
                        filter1[i][j][k][l] = (i+j+k+l) % 5;
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
            bias1[i] = i%3;
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
                        filter2[i][j][k][l] = (i+j+k+l) % 5;
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
            bias2[i] = i%3;
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
        for(int j = 0; j < H2; j+=4)
        {
            if(DEBUG)cout<<"output[][]["<<j<<"]\n";
            for(int k = 0; k < W2; k+=4)
            {
                for(int i = 0; i < D2; i++)
                {
                    int conH, conW, v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15;
                    //int w2, int h2, int d2
                    v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = v9 = v10 = v11 = v12 = v13 = v14 = v15 = bias1[i];
                    for(int iF = 0; iF < D1; iF++)
                    {
                        for(int jF = 0; jF <F; jF++)
                        {
                            for(int kF = 0; kF < F; kF++)
                            {
                                v0 += input[iF][jF+S*(j)][kF+S*(k)] * filter1[i][iF][jF][kF];
                                v4 += input[iF][jF+S*(j)][kF+S*(k+1)] * filter1[i][iF][jF][kF];
                                v8 += input[iF][jF+S*(j)][kF+S*(k+2)] * filter1[i][iF][jF][kF];
                                v12 += input[iF][jF+S*(j)][kF+S*(k+3)] * filter1[i][iF][jF][kF];
                                v1 += input[iF][jF+S*(j+1)][kF+S*(k)] * filter1[i][iF][jF][kF];
                                v5 += input[iF][jF+S*(j+1)][kF+S*(k+1)] * filter1[i][iF][jF][kF];
                                v9 += input[iF][jF+S*(j+1)][kF+S*(k+2)] * filter1[i][iF][jF][kF];
                                v13 += input[iF][jF+S*(j+1)][kF+S*(k+3)] * filter1[i][iF][jF][kF];
                                v2 += input[iF][jF+S*(j+2)][kF+S*(k)] * filter1[i][iF][jF][kF];
                                v6 += input[iF][jF+S*(j+2)][kF+S*(k+1)] * filter1[i][iF][jF][kF];
                                v10 += input[iF][jF+S*(j+2)][kF+S*(k+2)] * filter1[i][iF][jF][kF];
                                v14 += input[iF][jF+S*(j+2)][kF+S*(k+3)] * filter1[i][iF][jF][kF];
                                v3 += input[iF][jF+S*(j+3)][kF+S*(k)] * filter1[i][iF][jF][kF];
                                v7 += input[iF][jF+S*(j+3)][kF+S*(k+1)] * filter1[i][iF][jF][kF];
                                v11 += input[iF][jF+S*(j+3)][kF+S*(k+2)] * filter1[i][iF][jF][kF];
                                v15 += input[iF][jF+S*(j+3)][kF+S*(k+3)] * filter1[i][iF][jF][kF];
                            }
                        }
                    }
                    /*v0=sum1(k,j,i);
                    v1=sum1(k,j+1,i);
                    v2=sum1(k,j+2,i);
                    v3=sum1(k,j+3,i);
                    v4=sum1(k+1,j,i);
                    v5=sum1(k+1,j+1,i);
                    v6=sum1(k+1,j+2,i);
                    v7=sum1(k+1,j+3,i);
                    v8=sum1(k+2,j,i);
                    v9=sum1(k+2,j+1,i);
                    v10=sum1(k+2,j+2,i);
                    v11=sum1(k+2,j+3,i);
                    v12=sum1(k+3,j,i);
                    v13=sum1(k+3,j+1,i);
                    v14=sum1(k+3,j+2,i);
                    v15=sum1(k+3,j+3,i);*/
                    for(int ii = 0; ii < D3; ii++)
                    {
                        for(int jj = -F+1; jj <1; jj+=S)
                        {
                            for(int kk = -F+1; kk <1; kk+=S)
                            {
                                conH = jj+j+P;
                                conW = kk+k+P;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v0*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v4*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v8*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v12*filter2[ii][i][-jj][-kk];
                                }
                                conH+=1;
                                conW-=3;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v1*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v5*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v9*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v13*filter2[ii][i][-jj][-kk];
                                }
                                conH+=1;
                                conW-=3;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v2*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v6*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v10*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v14*filter2[ii][i][-jj][-kk];
                                }
                                conH+=1;
                                conW-=3;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v3*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v7*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v11*filter2[ii][i][-jj][-kk];
                                }
                                conW+=1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=v15*filter2[ii][i][-jj][-kk];
                                }
                            }
                        }
                    }
                    /*//cout<<"Computing, "<<i<<","<<j+P<<","<<k+P<<endl;
                    //cout<<"Accessing,"<<endl;
                    int val1 = sum1(k,j,i),val2 = sum1(k+1,j,i),val3 = sum1(k+1,j+1,i),val4 = sum1(k,j+1,i), conH, conW;
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
                                    output2[ii][conH][conW]+=val1*filter2[ii][i][-jj][-kk];
                                    //cout<<"Layer2,"<<ii<<","<<conH<<","<<conW<<endl;
                                    //cout<<"Filter2,"<<ii<<","<<i<<","<<-jj<<","<<-kk<<endl;
                                }
                                conW = conW+1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=val2*filter2[ii][i][-jj][-kk];
                                }
                                conH = conH+1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=val3*filter2[ii][i][-jj][-kk];
                                }
                                conW = conW-1;
                                if(P<=conW && conW<W3+P && P<=conH && conH<H3+P){
                                    output2[ii][conH][conW]+=val4*filter2[ii][i][-jj][-kk];
                                }
                            }
                        }
                    }*/
                }
                if(DEBUG)cout<<"\n";
            }
        }
        for(int i = 0; i < D3; i++)
        {
            for(int j = 0; j < H3; j++)
            {
                for(int k = 0; k < W3; k++)
                {
                    output2[i][j+P][k+P] += bias2[i];
                }
            }
        }
    }
    bool checkOutput(int ***output)
    {
        for(int i = 0; i < D3; i++)
        {
            for(int j = 0; j < H3; j++)
            {
                for(int k = 0; k < W3; k++)
                {
                    cout<<output2[i][j+P][k+P]<<","<<output[i][j+P][k+P]<<endl;
                    if(output2[i][j+P][k+P]!=output[i][j+P][k+P])return false;
                }
            }
        }
        return true;
    }
    bool checkInput(int ***input)
    {
        for(int i = 0; i < D1; i++)
        {
            for(int j = 0; j < H1+2*P; j++)
            {
                for(int k = 0; k < W1+2*P; k++)
                {
                    if(this->input[i][j][k] != input[i][j][k])return false;
                }
            }
        }
        return true;
    }
    bool checkFilter(int ****filter1, int ****filter2)
    {
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < D1; j++)
            {
                for(int k = 0; k < F; k++)
                {
                    for(int l = 0; l < F; l++)
                    {
                        if(filter1[i][j][k][l] != this->filter1[i][j][k][l])return false;
                    }
                }
            }
        }
        for(int i = 0; i < K; i++)
        {
            for(int j = 0; j < D2; j++)
            {
                for(int k = 0; k < F; k++)
                {
                    for(int l = 0; l < F; l++)
                    {
                        if(filter2[i][j][k][l] != this->filter2[i][j][k][l])return false;
                    }
                }
            }
        }
        return true;
    }
    bool checkBias(int *bias1, int *bias2)
    {
        for(int i = 0; i < K; i++)
        {
            if(bias1[i]!=this->bias1[i])return false;
        }
        for(int i = 0; i < K; i++)
        {
            if(bias2[i]!=this->bias2[i])return false;
        }
        return true;
    }
};

class Convolution
{
private:
    int ***input;
    int ****filter;
    int *bias;
    int ***output;
    int W1, H1, D1;
    int K, F, S, P;
    int W2, H2, D2;
public:
    Convolution(int W1, int H1, int D1, int K, int F, int S, int P, int ***outputFromPrevLayer=NULL)
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
        }
        if(outputFromPrevLayer == NULL)
        {
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
                            input[i][j][k] = (i+j+k) % 10;
                        }
                        else {
                            input[i][j][k] = 0;
                        }
                    }
                }
            }
        } else {
            this->input = outputFromPrevLayer;
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
        output = new int **[D2];
        for(int i = 0; i < D2; i++)
        {
            output[i] = new int *[H2+2*P];
            for(int j = 0; j < H2 + 2*P; j++)
            {
                output[i][j] = new int[W2+2*P];
                for(int k = 0; k < W2+2*P; k++)
                {
                    output[i][j][k] = 0;
                }
            }
        }
        filter = new int ***[K];
        for(int i = 0; i < K; i++)
        {
            filter[i] = new int **[D1];
            for(int j = 0; j < D1; j++)
            {
                filter[i][j] = new int *[F];
                for(int k = 0; k < F; k++)
                {
                    filter[i][j][k] = new int[F];
                    for(int l = 0; l < F; l++)
                    {
                        filter[i][j][k][l] = (i+j+k+l) % 5;
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
                            cout<<filter[i][j][k][l]<<" ";
                        }
                        cout<<endl;
                    }
                }
            }
        }
        bias = new int[K];
        for(int i = 0; i < K; i++)
        {
            bias[i] = i%3;
        }
        if(DEBUG)
        {
            cout<<"Bias"<<endl;
            for(int i = 0; i < K; i++)
            {
                cout << bias[i]<<" ";
            }
            cout<<endl;
        }
    }
    int sum(int w2, int h2, int d2)
    {
        int retValue = bias[d2];
        for(int j = 0; j <F; j++)
        {
            for(int k = 0; k < F; k++)
            {
                for(int i = 0; i < D1; i++)
                {
                    if(DEBUG)cout<<i<<j+S*h2<<k+S*w2<<i<<j<<k<<" ";
                    //cout<<i<<","<<j+S*h2<<","<<k+S*w2<<endl;
                    retValue += input[i][j+S*h2][k+S*w2] * filter[d2][i][j][k];
                }
            }
        }
        if(DEBUG)cout<<endl;
        return retValue;
    }
    void run()
    {
        //cout<<"start convolution\n";
        for(int i = 0; i < D2; i++)
        {
            if(DEBUG)cout<<"output[][]["<<i<<"]\n";
            for(int j = 0; j < H2; j++)
            {
                for(int k = 0; k < W2; k++)
                {
                    //cout<<"Computing, "<<i<<","<<j+P<<","<<k+P<<endl;
                    //cout<<"Accessing,"<<endl;
                    output[i][j+P][k+P] = sum(k,j,i);
                    if(DEBUG)cout<<output[i][j+P][k+P]<<" ";
                }
                if(DEBUG)cout<<"\n";
            }
        }
    }
    int *** getOutput(int *W2, int *H2, int *D2)
    {
        *W2 = this->W2;
        *H2 = this->H2;
        *D2 = this->D2;
        return output;
    }
    void copyBack(int ***outputFromPrevLayer)
    {
        for(int i = 0; i < D1; i++)
        {
            for(int j = 0; j < H1+2*P; j++)
            {
                for(int k = 0; k <W1+2*P; k++)
                {
                    /*if(P<=k && k<W1-1+2*P && P<=j && j<H1-1+2*P){
                        input[i][j][k] = rand() % 10;
                    }
                    else {
                        input[i][j][k] = 0;
                    }*/
                    input[i][j][k] = outputFromPrevLayer[i][j][k];
                }
            }
        }
    }
    int ***getInput()
    {
        return this->input;
    }
    int ****getFilter()
    {
        return this->filter;
    }
    int *getBias()
    {
        return this->bias;
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
    Convolution firstLayer(W1, H1, D1, K, F, S, P);
    int ***firstLayerOutput = firstLayer.getOutput(&W1, &H1, &D1);
    Convolution secondLayer(W1, H1, D1, K, F, S, P, firstLayerOutput);
    struct timeval tv;
    gettimeofday(&tv,NULL);
    long start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    fusion.run();
    gettimeofday(&tv,NULL);
    long end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    long elapsed = end - start;
    cout<<"s, "<<(1.0 * elapsed)/1000000<<endl;
    /*gettimeofday(&tv,NULL);
    start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    firstLayer.run();
    secondLayer.run();
    gettimeofday(&tv,NULL);
    end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    elapsed = end - start;
    cout<<"s, "<<(1.0 * elapsed)/1000000<<endl;
    firstLayerOutput = secondLayer.getOutput(&W1, &H1, &D1);
    cout<<fusion.checkInput(firstLayer.getInput())<<endl;
    cout<<fusion.checkFilter(firstLayer.getFilter(), secondLayer.getFilter())<<endl;
    cout<<fusion.checkBias(firstLayer.getBias(), secondLayer.getBias())<<endl;
    cout<<fusion.checkOutput(firstLayerOutput)<<endl;*/
    return 0;
}

