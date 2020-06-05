#include<iostream>
#include<algorithm>
#include "QR_factorization.h"

using namespace std;

vector<double> QRfactor::work(double* A) const{
    int n = a, m = x.size();
    double G[(n+1)*(n+1)], c[n+1];
    vector<double> anwser;
    Generatematrix(G, c, m, n+1);
    for(auto i = 0; i < ((n+1)*(n+1)); ++i){
        A[i] = G[i];
    }
    Equationsolver(G, c, n+1);
    for_each(c, c + n+1, [&anwser](double i){ anwser.push_back(i);});
    return anwser;
}

void QRfactor::Generatematrix(double* G, double* c, const int m, const int n) const{
    double** A = new double*[m],
        **Q = new double*[m];
    for(auto i = 0; i < m; i++){
        A[i] = new double[n];
        Q[i] = new double[m];
    }
    for(auto i = 0; i < m; ++i){
        for(auto j = 0; j < n; ++j){
            A[i][j] = pow(x[i], j);
        }
    }
    Maqr(A, Q, m, n);
    for(auto i = 0; i < n; ++i){
        for(auto j = 0; j < n; ++j){
            G[j*n + i] = A[i][j];
        }
    }
    for(auto i = 0; i < n; ++i){
        double sum = 0;
        for(auto j = 0; j < m; ++j){
            sum += Q[j][i] * y[j];
        }
        c[i] = sum;
    }
}

void QRfactor::Maqr(double** A, double** Q, const int m, const int n) const{
    int i, j, k, nn, jj;
	double u, alpha, w, t;
	//double** Q = new double*[m];   //动态分配内存空间
	//for (i = 0; i<m; i++) Q[i] = new double[m];
	if (m < n)
	{
		cout << "\nQR分解失败！" << endl;
		exit(1);
	}
	//保证列数<=行数
	for (i = 0; i <= m - 1; i++)
		for (j = 0; j <= m - 1; j++)
		{
			Q[i][j] = 0.0;
			if (i == j) Q[i][j] = 1.0;
		}
    //初始的Q矩阵就是一个单位的m阶方阵
    nn = n;
	if (m == n) nn = m - 1;
	for (k = 0; k <= nn - 1; k++)//在大循环k：0~m当中，进行H矩阵的求解，左乘Q，以及左乘A
	{

		u = 0.0;
		for (i = k; i <= m - 1; i++)
		{
			w = fabs(A[i][k]);
			if (w > u) u = w;
		}
		alpha = 0.0;
		for (i = k; i <= m - 1; i++)
		{
			t = A[i][k] / u; alpha = alpha + t * t;
		}
		if (A[k][k] > 0.0) u = -u;
		alpha = u * sqrt(alpha);
		if (fabs(alpha) + 1.0 == 1.0)
		{
			cout << "\nQR分解失败！" << endl;
			exit(1);
		}

		u = sqrt(2.0*alpha*(alpha - A[k][k]));
		if ((u + 1.0) != 1.0)
		{
			A[k][k] = (A[k][k] - alpha) / u;
			for (i = k + 1; i <= m - 1; i++) A[i][k] = A[i][k] / u;
			
			//以上就是H矩阵的求得，实际上程序并没有设置任何数据结构来存储H矩
			//阵，而是直接将u向量的元素赋值给原A矩阵的原列向量相应的位置，这样做
			//这样做是为了计算左乘矩阵Q和A
			for (j = 0; j <= m - 1; j++)
			{
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A[jj][k] * Q[jj][j];
				for (i = k; i <= m - 1; i++)
					Q[i][j] = Q[i][j] - 2.0*t*A[i][k];
			}
//左乘矩阵Q，循环结束后得到一个矩阵，再将这个矩阵转置一下就得到QR分解中的Q矩阵
//也就是正交矩阵

			for (j = k + 1; j <= n - 1; j++)
			{
				t = 0.0;
				for (jj = k; jj <= m - 1; jj++)
					t = t + A[jj][k] * A[jj][j];
				for (i = k; i <= m - 1; i++)
					A[i][j] = A[i][j] - 2.0*t*A[i][k];
			}
			//H矩阵左乘A矩阵，循环完成之后，其上三角部分的数据就是上三角矩阵R
			A[k][k] = alpha;
			for (i = k + 1; i <= m - 1; i++)  A[i][k] = 0.0;
		}
	}
	for (i = 0; i <= m - 2; i++)
		for (j = i + 1; j <= m - 1; j++)
		{
			t = Q[i][j]; Q[i][j] = Q[j][i]; Q[j][i] = t;
		}
	//QR分解完毕，然后在函数体里面直接将Q、R矩阵打印出来
    /*
	for (i = 0; i<m; i++)
	{
		for (j = 0; j<n; j++)
		{
			cout << "    " << Q[i][j];
		}
		cout << endl;
	}
	cout << endl;
	for (i = 0; i<m; i++)
	{
		for (j = 0; j<n; j++)
		{
			cout << "    " << A[i][j];
		}
		cout << endl;
	}
    */
}
