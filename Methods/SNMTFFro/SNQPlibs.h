#ifndef __LIBS__
//#include <iostream>
#include <cmath>
#include <cstring>
#include <algorithm>
#include "pthread.h"
#include "math.h"
#include "mex.h"
#include <time.h>
#include <cmath>
#include <stdlib.h>   
#define MAX_SUPPORTED_THREADS 100

#define For(i,n) for(int i=0; i<n; i++)
//#define max(a, b) (a > b ?  a : b)
//#define abs(x) (x >= 0 ? x : -x)
#define eps 1e-10

using namespace std;

//inner product of two vectors
double dotVVFrom(double *a, double *b, int from, int n){
    double res = 0; 
    double *ai, *bi, *end=a+n;
    for (ai=(a+from), bi=(b+from); ai!=end; ai++,bi++) res += (*ai)*(*bi);
    return res;
}

//assign vector x = vector y
void assign(double *x, double *y, int k){
    For(i, k) x[i] = y[i];
}

//inner product of two vectors
double dotVV(double *a, double *b, int n){
    return dotVVFrom(a, b, 0, n);
}

//sum of two vectors
void addVV(double *a, double *b, int n, double *res){
    double *end = res+n;
    for(double *ir=res, *ia=a, *ib=b; ir!=end; ir++, ia++, ib++)
        *ir = *ia + *ib;
}

void addkV(double alpha, double *a, int n, double *res){
    double *end=a+n;
    for(double* it=a, *ires=res; it != end; it++,ires++) (*ires) += alpha*(*it);
}

//subtract of two vectors
void subVV(double *a, double *b, int n, double *res){
    For(i,n) res[i] = a[i] - b[i];
}

//a product of matrix and vector
void dotMV(double **a, double *b, int m, int n, double *res){
    For(i,m) res[i] = dotVV(a[i], b, n);
}

//a product of matrix and vector
void dotVM(double **a, double *b, int m, int n, double *res){
    For(i,m) res[i] = dotVV(a[i], b, n);
}

//a product of v^TQv
double dotVMV(double *x, double **M, int k){
    double res = 0, sub;
    for (int i = 0; i < k; i++) {
        res += x[i] * x[i] * M[i][i];
        sub = 0;
        sub = dotVVFrom(M[i]+i+1, x+i+1, i+1, k-i-1);
        //for (int j = i + 1; j < k; j++)
        //    sub += M[i][j] * x[j];
        res += 2 * x[i] * sub;
    }
    return res;
}

//Remove negative elements of a vector
void removeVVNegative(double* a, int n, double *b){
    For(i, n) b[i] = (a[i] >= 0 ? a[i] : 0);
}

//Remove KKT-satisfied elements of a vector
void removeKTTElements(double* df, double* x, int n, double *dbarf){
    For(i, n) dbarf[i] = (df[i] < 0 || x[i] > 0) ? df[i] : 0;
}


//fro norm ||a||^2_2 
double square(double *a, int n){
    double res = 0;
    For(i, n) res += a[i]*a[i];
    return res;
}

//get index of max element
int maxV(double *a, int n){
    int k = 0;
    For(i, n) k = (a[i] > a[k] ? i : k);
    return k;
}

//convert double* into double**
double **convert(double *a, int m, int n){
    double **res = (double **)malloc(sizeof(double*)*m);
    For(i,m) res[i] = &a[i*n];
    return res;
}

//create a vector
double *createV(int n){
    double *res = (double *)malloc(sizeof(double)*n);
    For(i,n) res[i] = 0;
    return res;
}

//scale a vector
void scaleV(double *a, int n, double s){
    double* end=a+n;
    for(double* iter=a; iter!=end; iter++) 
        (*iter) *= s;
}

double getMaxReduce(double **H, double *h, int k, double* df, double* x){
    double dx, maxValue=-1e30, dfx;
    int pos = -1, v;
    For(i, k)
        if ((df[i] < 0 || x[i] > 0) && H[i][i] > 0){
            dx = max(0.0, x[i] - df[i]) - x[i];
            dfx = (-0.5 * dx - df[i]) * dx;
            if (maxValue < dfx){
                pos = i;
                maxValue = dfx;
            }
        }
    return maxValue;
}

double sum(double *x, int k){
    double res = 0;
    For(i, k) res += x[i];
    return res;
}

int getPosMin(double *x, int k){
    double vMin = x[0];
    int pos = 0;
    For(i, k)
        if (x[i] < vMin){
            vMin = x[i];
            pos = i;
        }
    return pos;
}

class thread_data {
public:
    double **HH;
    double **hh;
    int k;
    double tolorance;
    int maxIter;
    int verbose;
    double **xx;
    double *iter;
    int start;
    int end;
    double *sqrtH;
    double *iSqrtH;
    
    double *df;
    double *dbarf;
    double *newX;
    double *saveDf;
    double *saveX;
    double *dx;
    
    thread_data(){
        df = createV(k);
        dbarf = createV(k);
        newX = createV(k);

        saveDf = createV(k);
        saveX = createV(k);
        dx = createV(k);
    }
    
    thread_data(double **HH, double **hh, int k, 
            double tolorance, int maxIter, int verbose,
            double **xx, double *iter, int start, int end, double *sqrtH, double *iSqrtH){
        this->HH = HH;
        this->hh = hh;
        this->k = k;
        this->tolorance = tolorance;
        this->maxIter = maxIter;
        this->verbose = verbose;
        this->xx = xx;
        this->iter = iter;
        this->start = start;
        this->end = end;
        this->sqrtH = sqrtH;
        this->iSqrtH = iSqrtH;
        
        df = createV(k);
        dbarf = createV(k);
        newX = createV(k);

        saveDf = createV(k);
        saveX = createV(k);
        dx = createV(k);
    }
    
    ~thread_data(){
        free(df);
        free(dbarf);
        free(newX);

        free(saveDf);
        free(saveX);
        free(dx);
    }
};

void solveSNQP(double** Q, double *q, int k, double tolorance, int maxIter, int verbose,
        double *x, double* iter, double *stopMax, double* sqrtH, double* iSqrtH, thread_data* args){
    double *df = args->df;
    double *Qx = args->dbarf;
    //double *qx = args->newX;
    
    double *saveDf = args->saveDf;
    double *saveX = args->saveX;
    double *dx = args->dx;
    
    double s = sum(x, k), alpha, qx = 0, init;
    
    if (abs(s) < 0.001) {
        For(i, k) x[i] = 0;
        int pos = 0; 
        double vmin = Q[0][0]/2.0 + q[0];
        For(i, k)
            if (vmin > Q[i][i]/2 + q[i])
                pos = i, vmin = Q[i][i]/2 + q[i];
        For(i, k) Qx[i] = df[i] = x[i] = 0;
        x[pos] = 1.0;
        init = Q[pos][pos]/2 + q[pos];
    } /*else {
       dotMV(Q, x, k, k, Qx);
       addVV(Qx, q, k, df);
       qx = dotVV(q, x, k);
       init = dotVV(x, Qx, k)/2 + qx; 
       
       For(i, k) x[i] = 0;
        int pos = 0; 
        double vmin = Q[0][0]/2.0 + q[0];
        For(i, k)
            if (vmin > Q[i][i]/2 + q[i])
                pos = i, vmin = Q[i][i]/2 + q[i];
        For(i, k) Qx[i] = df[i] = x[i] = 0;
        x[pos] = 1.0;
    } /**/
    dotMV(Q, x, k, k, Qx);
    addVV(Qx, q, k, df);
    qx = dotVV(q, x, k);
    //addVV(Qx, q, k, df);
    
    double den, a, error, xQx;
    
    *iter = 0;
    if (abs (sum(x, k) - 1.0) > 1e-5)
        mexPrintf(">> sum(x) = %f,   %f \n", sum(x, k), s);
   
    For(it, 1){
        *iter = *iter + 1;
        //assign(saveDf, df, k);
        //assign(saveX, x, k);
        double vmax, dfa;
        int maxIter = min(4*k, 100);
        for(int t=0; t < maxIter; t++){
            int v = 0, pos; 
            double xdf = dotVV(df, x, k);
            xQx = dotVV(x, Qx, k);
            vmax = 0;
            For(i, k) {
                //den = Q[i][i] - 2*Qx[i] + xQx;
                dfa = (Qx[i] - xQx + q[i] - qx); 
                //if (alpha != alpha || abs(x[i] - 1.0) < eps) continue;
                if ((dfa < 0 || (dfa > 0 && x[i] > 0)) && abs(dfa) > vmax) 
                    vmax = abs(dfa), v = i;
            }
            pos = abs(v);
            
            if (abs(x[pos] - 1.0) < eps)
                break;
            
            den = Q[pos][pos] - 2*Qx[pos] + xQx;

            alpha = -(Qx[pos] - xQx + q[pos] - qx) / den;

            if (alpha != alpha || abs(alpha) < eps) break;
            alpha = min(alpha, 1.0);
            alpha = max(alpha, -x[pos]/(1.0 - x[pos]));
            
            if (t == 0)
                init = abs(vmax);
            else {
                if (abs(vmax) < *stopMax * 0.1) {
                   *stopMax = max(*stopMax, vmax);
                   break;
                }
                if ( abs(vmax) < init * 0.1) 
                    break;
                init = max(init, abs(vmax));
            }
            
            scaleV(Qx, k, 1 - alpha);
            scaleV(x, k, 1 - alpha);
            
            qx = (1 - alpha) * qx + alpha * q[pos];
         
            x[pos] = max(0.0, x[pos] + alpha);
            addkV(alpha, Q[pos], k, Qx);
            
            addVV(Qx, q, k, df);
            //if ((dotVV(df, x, k) + dotVV(x, q, k))/2.0 < init*0.01) break; 
            //if (vmin > 0)
            //if (abs (sum(x, k) - 1.0) > 1e-5 || alpha < 0 || v < 0)
            //mexPrintf("%d sum(x)=%f p=%d, a=%f : %f %f %f\n", t, sum(x, k), v, alpha, 
            //        dotVV(x, Qx, k)/2.0 + qx, 
            //        (dotVV(df, x, k) + dotVV(x, q, k))/2.0, 
            //         dotVMV(x, Q, k)/2.0 + dotVV(x, q, k));
        }
        *stopMax = max(*stopMax, vmax);
    }
}


void* thread_iter_SNQP(void* thread_args){
    thread_data *args = (thread_data *) thread_args;
    double stopMax = 0;
    for (int i=args->start; i < args->end; i++) {
        double it = 0;
        solveSNQP
        (args->HH, args->hh[i], args->k, args->tolorance, args->maxIter, args->verbose, args->xx[i], &it, &stopMax, args->sqrtH, args->iSqrtH, args);
        *(args->iter) += it;
    }
    return NULL;
}


void print(double* a, int k){
    For(i,k) mexPrintf("%.4f ", a[i]);
    mexPrintf("\n");
}

void usage()
{
	printf("Error calling doiter.\n");
	printf("Usage: Wnew = doiter(GW, HH^T, W, tol, maxinner)\n");
}

void getArray(const mxArray* a, double **res, int *m, int *n){
    *res = mxGetPr(a);
    *m = mxGetM(a);
    *n = mxGetN(a);
}

double getDouble(const mxArray* a){
    double *values = mxGetPr(a);
    return values[0];
}

#endif
