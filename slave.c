#include <slave.h>
#include "matrix_def.h"
#include <stdio.h>
// #include <assert.h>

#define CORE_NUM       64
#define N              200000
#define WRITE_BUF_SIZE 1024
#define READ_BUF_SIZE  1024

#define isnan(x) (!((x) == (x)))
#define min(a, b) ((a) < (b) ? (a) : (b))

typedef struct{
    struct CsrMatrix csr_matrix;
    double *x;
    double *b;
} Para;

__thread_local crts_rply_t dma_rply[2] = {0, 0}; // 局存变量，dma传输回答字地址
__thread_local crts_rply_t dma_rply3 = 0; // 局存变量，dma传输回答字地址
__thread_local unsigned int D_COUNT[2] = {0, 0}; // 局存变量，需等待完成的非阻塞DMA次数
__thread_local unsigned int D_COUNT3 = 0; // 局存变量，需等待完成的非阻塞DMA次数
__thread_local int row_off[N / CORE_NUM * 2];
// __thread_local int *row_off;
__thread_local int cols[2 * READ_BUF_SIZE];
__thread_local double data[2 * READ_BUF_SIZE];
__thread_local double write_buf[2 * WRITE_BUF_SIZE];

__thread_local_share double X[N];
// __thread_local_share double B[N];

int Find(int l, int n, int offset) {
    // assert(row_off[l + 1 - offset] - row_off[l - offset] <= READ_BUF_SIZE);
    // if (row_off[l + 1] - row_off[l] > READ_BUF_SIZE) ;
    int r = l;
    while (r < n && row_off[r + 1 - offset] - row_off[l - offset] <= READ_BUF_SIZE) ++r;
    return r;
}

void slave_spmv(void *arg) {
    Para *para = arg;
    int id = CRTS_tid; // 从核id
    const struct CsrMatrix csr_matrix = para->csr_matrix;
    double *x = para->x;
    double *b = para->b;
    const int n = csr_matrix.rows;

    athread_memcpy_sldm(X, x, n * sizeof(double), MEM_TO_LDM);

    int block = (csr_matrix.rows + CORE_NUM - 1) / CORE_NUM;
    int l = block * id, r = min(block * (id + 1), csr_matrix.rows);
    int wboffset = 0, size = 0;
    int bs = l;
    // assert((r - l + 1) <= N / CORE_NUM * 2);
    // row_off = CRTS_pldm_malloc((r - l + 1) * sizeof(int));
    // if (row_off == NULL) printf("%d\n", (r - l + 1) * sizeof(int) / 1024);
    // return;
    // assert(row_off != NULL);
    CRTS_dma_iget(row_off, csr_matrix.row_off + l, (r - l + 1) * sizeof(int), &dma_rply3);
    D_COUNT3++;
    CRTS_dma_wait_value(&dma_rply3, D_COUNT3);
    int cur = 0;
    int ll[2] = {0, 0}, rr[2] = {0, 0};
    ll[cur] = l; rr[cur] = Find(l, r, l);
    // assert(rr[cur] == r);
    // assert(ll[cur] != rr[cur]);
    CRTS_dma_iget(cols + cur * READ_BUF_SIZE, csr_matrix.cols + row_off[ll[cur] - l], (row_off[rr[cur] - l] - row_off[ll[cur] - l]) * sizeof(int), &dma_rply[cur]);
    CRTS_dma_iget(data + cur * READ_BUF_SIZE, csr_matrix.data + row_off[ll[cur] - l], (row_off[rr[cur] - l] - row_off[ll[cur] - l]) * sizeof(double), &dma_rply[cur]);
    D_COUNT[cur] += 2;
    ll[cur ^ 1] = rr[cur]; rr[cur ^ 1] = Find(ll[cur ^ 1], r, l);
    if (ll[cur ^ 1] < rr[cur ^ 1]) {
        CRTS_dma_iget(cols + (cur ^ 1) * READ_BUF_SIZE, csr_matrix.cols + row_off[ll[cur ^ 1] - l], (row_off[rr[cur ^ 1] - l] - row_off[ll[cur ^ 1] - l]) * sizeof(int), &dma_rply[cur ^ 1]);
        CRTS_dma_iget(data + (cur ^ 1) * READ_BUF_SIZE, csr_matrix.data + row_off[ll[cur ^ 1] - l], (row_off[rr[cur ^ 1] - l] - row_off[ll[cur ^ 1] - l]) * sizeof(double), &dma_rply[cur ^ 1]);
        D_COUNT[cur ^ 1] += 2;
    }
    CRTS_dma_wait_value(&dma_rply[cur], D_COUNT[cur]);
    // for (int i = id; i < csr_matrix.rows; i += CORE_NUM) {
    for (int i = l; i < r; ++i) {
        int start = row_off[i - l] - row_off[ll[cur] - l];
        int num = row_off[i + 1 - l] - row_off[i - l];
        double result = .0;
        for(int j = 0; j < num; j++) {
            int k = cols[cur * READ_BUF_SIZE + start + j];
            // assert(k == csr_matrix.cols[csr_matrix.row_off[i]+j]);
            // assert(data[cur * READ_BUF_SIZE + start + j] == csr_matrix.data[csr_matrix.row_off[i]+j]);
            result += (k < n ? X[k] : x[k]) * data[cur * READ_BUF_SIZE + start + j];
            // int k = csr_matrix.cols[start+j];
            // result += (k < n ? X[k] : x[k]) * csr_matrix.data[start+j];
        }
        if (i + 1 == rr[cur]) {
            ll[cur] = rr[cur ^ 1]; rr[cur] = Find(ll[cur], r, l);
            if (ll[cur] < rr[cur]) {
                CRTS_dma_iget(cols + cur * READ_BUF_SIZE, csr_matrix.cols + row_off[ll[cur] - l], (row_off[rr[cur] - l] - row_off[ll[cur] - l]) * sizeof(int), &dma_rply[cur]);
                CRTS_dma_iget(data + cur * READ_BUF_SIZE, csr_matrix.data + row_off[ll[cur] - l], (row_off[rr[cur] - l] - row_off[ll[cur] - l]) * sizeof(double), &dma_rply[cur]);
                D_COUNT[cur] += 2;
            }
            cur ^= 1;
            CRTS_dma_wait_value(&dma_rply[cur], D_COUNT[cur]);
        }
        // B[i] = result;
        // b[i] = result;
        write_buf[wboffset + size++] = result;
        if (size == WRITE_BUF_SIZE) {
            CRTS_dma_wait_value(&dma_rply3, D_COUNT3);
            CRTS_dma_iput(b + bs, write_buf + wboffset, WRITE_BUF_SIZE * sizeof(double), &dma_rply3);
            D_COUNT3++;
            wboffset ^= WRITE_BUF_SIZE;
            bs = i + 1;
            size = 0;
        }
    }
    // CRTS_ssync_sldm();
    // athread_memcpy_sldm(b, B, n * sizeof(double), LDM_TO_MEM);
    if (size > 0) {
        CRTS_dma_iput(b + bs, write_buf + wboffset, size * sizeof(double), &dma_rply3);
        D_COUNT3++;
    }
    // CRTS_pldm_free(row_off, (r - l + 1) * sizeof(int));
    CRTS_dma_wait_value(&dma_rply3, D_COUNT3);
}