#include <athread.h>
#include "matrix_def.h"
#include "host.h"

extern "C" {
    void slave_spmv(void *arg);
}
typedef struct{
    CsrMatrix csr_matrix;
    double *x;
    double *b;
} Para;
Para para;

void spmv(const CsrMatrix &csr_matrix, double *x, double *b) {
    // 实现样例
    // for(int i = 0; i < csr_matrix.rows; i++) {
    //     int start = csr_matrix.row_off[i];
    //     int num = csr_matrix.row_off[i+1] - csr_matrix.row_off[i];
    //     double result = .0;
    //     for(int j = 0; j < num; j++) {
    //         result += x[csr_matrix.cols[start+j]] * csr_matrix.data[start+j];
    //     }
    //     b[i]=result;
    // }
    para.csr_matrix = csr_matrix;
    para.x = x;
    para.b = b;
    CRTS_init();
    athread_spawn(slave_spmv, &para);
    athread_join();
    athread_halt();
    // int cnt = 0;
    // for(int i = 0; i < csr_matrix.rows; i++)
    //     if (b[i] == 0) ++cnt;
    // printf("%d/%d\n", cnt, csr_matrix.rows);
}
