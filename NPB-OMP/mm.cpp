#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mm_rdma.cpp"



#define MATRIX_SIZE 400

// Assuming matrices are square and of size MATRIX_SIZE
double A[MATRIX_SIZE][MATRIX_SIZE];
double B[MATRIX_SIZE][MATRIX_SIZE];
double C[MATRIX_SIZE][MATRIX_SIZE];


void matrix_multiply() {
    #pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            C[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main() {
    //struct rdma_cm_id *rdma_id; // This should be properly initialized for RDMA communication
    struct resources res;
	int rc = 1;
	char temp_char;
	resources_init(&res);
	if (resources_create(&res)) {
		fprintf(stderr, "failed to create resources\n");
	}

	if (connect_qp(&res)) {
		fprintf(stderr, "failed to connect QPs\n");
	}

	if (sock_sync_data(res.sock, 1, "W", &temp_char)) /* just send a dummy char back and forth */
	{
		fprintf(stderr, "sync error after RDMA ops\n");
		rc = 1;
	}

    
    uint64_t A = mymalloc(&res, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, 0);
    if (poll_completion(&res)) {
		fprintf(stderr, "poll completion failed\n");
	}

    uint64_t B = mymalloc(&res, sizeof(double) * MATRIX_SIZE * MATRIX_SIZE, 1600);
    if (poll_completion(&res)) {
		fprintf(stderr, "poll completion failed\n");
	}

    
    myread(&res, A, res.buf, MATRIX_SIZE * MATRIX_SIZE);
    if (poll_completion(&res)) {
		fprintf(stderr, "poll completion failed\n");
	}

    myread(&res, B, res.buf + 1600, MATRIX_SIZE * MATRIX_SIZE);
    if (poll_completion(&res)) {
		fprintf(stderr, "poll completion failed\n");
	}
	//srand(time(NULL));
	//fill_memory_with_random(res.buf, MATRIX_SIZE * MATRIX_SIZE);
	
	//fill_memory_with_random(res.buf + 16, MATRIX_SIZE * MATRIX_SIZE);

	for (int i = 0; i < MATRIX_SIZE; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			res.buf[i * MATRIX_SIZE + j] = i + j;
		}
		//printf("\n");
	}

	for (int i = 0; i < MATRIX_SIZE; i++) {
		for (int j = 0; j < MATRIX_SIZE; j++) {
			res.buf[1600 + i * MATRIX_SIZE + j] = i + j;
		}
		//printf("\n");
	}

	#pragma omp parallel for
    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            C[i][j] = 0;
            for (int k = 0; k < MATRIX_SIZE; k++) {
                //C[i][j]+=A[i][k] * B[k][j];
                C[i][j] += res.buf[i * MATRIX_SIZE + k] * res.buf[1600 +  k * MATRIX_SIZE + j];
            }
        }
    }

	FILE *fp;
	fp = fopen("output.txt", "a"); 
	if (fp == NULL) {
    	perror("Error opening file");
	}

    for (int i = 0; i < MATRIX_SIZE; i++) {
        for (int j = 0; j < MATRIX_SIZE; j++) {
            fprintf(fp, "%f ", C[i][j]);
        }
        fprintf(fp, "\n");
    } 
    fclose(fp);
    
    if (sock_sync_data(res.sock, 1, "R", &temp_char)) /* just send a dummy char back and forth */
	{
		fprintf(stderr, "sync error before RDMA ops\n");
		rc = 1;
	}
    return 0;
}


// void fill_memory_with_random(double *buffer, size_t size) {
// 	srand(time(NULL));
//     for (size_t i = 0; i < size; i++) {
//         buffer[i] = rand() % 256; // 生成一个0到255之间的随机数
//     }
// }