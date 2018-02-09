/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: takh
 *
 * Created on January 31, 2018, 7:48 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

#define DEBUG 0

using namespace std;

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
                            input[i][j][k] = rand() % 10;
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
                        filter[i][j][k][l] = rand() % 5;
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
            bias[i] = rand()%3;
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
};



/*
 * 
 */
int main(int argc, char** argv) {
    ifstream inFile("input.txt");
    int W1, H1, D1, K, F, S, P;
    inFile >> W1 >> H1 >> D1;//5 5 3
    inFile >> K >> F >> S >> P;//2 3 2 1
    Convolution firstLayer(W1, H1, D1, K, F, S, P);
    int ***firstLayerOutput = firstLayer.getOutput(&W1, &H1, &D1);
    Convolution secondLayer(W1, H1, D1, K, F, S, P, firstLayerOutput);
    
    struct timeval tv;
    gettimeofday(&tv,NULL);
    long start = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    firstLayer.run();
    secondLayer.run();
    gettimeofday(&tv,NULL);
    long end = tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
    long elapsed = end - start;
    cout<<"s, "<<(1.0 * elapsed)/1000000<<endl;
    return 0;
}

