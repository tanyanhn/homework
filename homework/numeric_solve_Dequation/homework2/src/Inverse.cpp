#include "Inverse.h"
#include <memory.h>


namespace integrator_ty{

int SCalcInverseMatrix(Real* pDst, const Real* pSrc, int dim)
{
    int nRetVal = 0;
    if (pSrc == pDst)
    {
        return -1;
    }

    int* ipiv = new int[dim * dim];
    Real* pSrcBak = new Real[dim * dim];  // LAPACKE_sgesv会覆盖A矩阵，因而将pSrc备份
    memcpy(pSrcBak, pSrc, sizeof(Real)* dim * dim);

    memset(pDst, 0, sizeof(Real)* dim * dim);
    for (int i = 0; i < dim; ++i)
    {
        // LAPACKE_sgesv函数计算AX=B，当B为单位矩阵时，X为inv(A)
        pDst[i*(dim + 1)] = 1.0f;
    }

    // 调用LAPACKE_sgesv后，会将inv(A)覆盖到X（即pDst）中
    nRetVal = LAPACKE_dgesv(LAPACK_ROW_MAJOR, dim, dim, pSrcBak, dim, ipiv, pDst, dim);

    delete[] ipiv;
    ipiv = nullptr;

    delete[] pSrcBak;
    pSrcBak = nullptr;

    return nRetVal;
}

void Inverse::Solverequations(Real* A, Real* b){
    int n = dim;
    int ipiv[dim];
    int nrhs = 1;
    int info;
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs,
                              A, n, ipiv,
                              b, nrhs);
}

void print_matrix(lapack_int m, lapack_int n, double* a, lapack_int lda )
{
    lapack_int i, j;
    printf( "\n a \n");
    for( i = 0; i < m; i++ )
    {
        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
        printf( "\n" );
    }
}


void Inverse::operator()(Real* A){
    int m = dim;
    int n = dim;
    int lda = dim;
    int ipiv[dim];
    Real* d = new Real[dim*dim];
    int info;
    //print_matrix(m,n,A,lda);	
    //info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,m,n,A,lda,ipiv);
    //print_matrix(m,n,A,lda);
 
    //info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,m,A,lda,ipiv);
    //print_matrix(m,n,A,lda);
    SCalcInverseMatrix(d, A, dim);
    for(auto k =0; k < dim*dim; ++k){
        A[k] = d[k];
    }
    delete[] d;
    return ;

}

}