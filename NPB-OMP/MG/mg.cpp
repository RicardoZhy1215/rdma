/*
MIT License

Copyright (c) 2021 Parallel Applications Modelling Group - GMAP 
	GMAP website: https://gmap.pucrs.br
	
	Pontifical Catholic University of Rio Grande do Sul (PUCRS)
	Av. Ipiranga, 6681, Porto Alegre - Brazil, 90619-900

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

------------------------------------------------------------------------------

The original NPB 3.4.1 version was written in Fortran and belongs to: 
	http://www.nas.nasa.gov/Software/NPB/

Authors of the Fortran code:
	E. Barszcz
	P. Frederickson
	A. Woo
	M. Yarrow
	H. Jin

------------------------------------------------------------------------------

The serial C++ version is a translation of the original NPB 3.4.1
Serial C++ version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-SER

Authors of the C++ code: 
	Dalvan Griebler <dalvangriebler@gmail.com>
	Gabriell Araujo <hexenoften@gmail.com>
 	Júnior Löff <loffjh@gmail.com>

------------------------------------------------------------------------------

The OpenMP version is a parallel implementation of the serial C++ version
OpenMP version: https://github.com/GMAP/NPB-CPP/tree/master/NPB-OMP

Authors of the OpenMP code:
	Júnior Löff <loffjh@gmail.com>
	
*/
#include <infiniband/verbs.h>
#include "omp.h"
#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <unordered_map>
#include <endian.h>
#include <byteswap.h>
#include <getopt.h>
#include <algorithm>
#include <malloc.h>
#include <atomic>
#include <sys/time.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <math.h>
#define MAX_POLL_CQ_TIMEOUT 200000

using std::min;

/* structure of test parameters */
struct config_t
{
	char *dev_name; /* IB device name */
	char *server_name;	/* server host name */
	u_int32_t tcp_port;   /* server TCP port */
	int ib_port;		  /* local IB port to work with */
	int gid_idx;		  /* gid index to use */
};
/* structure to exchange data which is needed to connect the QPs */
struct cm_con_data_t
{
	uint64_t addr;   /* Buffer address */
	uint32_t rkey;   /* Remote key */
	uint32_t qp_num; /* QP number */
	uint16_t lid;	/* LID of the IB port */
	uint8_t gid[16]; /* gid */
} __attribute__((packed));

/* structure of system resources */
struct resources
{
	struct ibv_device_attr
		device_attr;
	/* Device attributes */
	struct ibv_port_attr port_attr;	/* IB port attributes */
	struct cm_con_data_t remote_props; /* values to connect to remote side */
	struct ibv_context *ib_ctx;		   /* device handle */
	struct ibv_pd *pd;				   /* PD handle */
	struct ibv_cq *cq;				   /* CQ handle */
	struct ibv_qp *qp;				   /* QP handle */
	struct ibv_mr *mr;				   /* MR handle for buf */
	double *buf;						   /* memory buffer pointer, used for RDMA and send
ops */
	int sock;						   /* TCP socket file descriptor */
	uint64_t atomic_compare;  // CAS 操作的比较值
    uint64_t atomic_swap;     // CAS 操作的交换值
    uint64_t atomic_add;      // F&A 操作的加值
};
struct config_t config = {
	"mlx5_2",  /* dev_name */
	"10.140.82.141",  /* server_name */
	19870, /* tcp_port */
	1,	 /* ib_port */
	0 /* gid_idx */};

void shuffle(size_t *array, size_t n);
static void resources_init(struct resources *res);
static int resources_create(struct resources *res);
static int sock_connect(const char *servername, int port);
int sock_sync_data(int sock, int xfer_size, char *local_data, char *remote_data);
static int poll_completion(struct resources *res);
static int poll_completion_2(struct resources *res);
static int modify_qp_to_init(struct ibv_qp *qp);
static int modify_qp_to_rtr(struct ibv_qp *qp, uint32_t remote_qpn, uint16_t dlid, uint8_t *dgid);
static int modify_qp_to_rts(struct ibv_qp *qp);
static int connect_qp(struct resources *res);
static int resources_destroy(struct resources *res);
uint64_t mymalloc(struct resources *res, size_t mem_size, size_t offset);
int mywrite(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size);
int myread(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size);


const size_t block_num = 1000;
const size_t num_threads = 1;
const size_t cpupool_mem_size = 4800 * 1024 / 100 * 100;
const size_t block_size =  cpupool_mem_size / num_threads;

const size_t mempool_mem_size = 1024ULL * 1024  * 1024 * 10;
           

#if __BYTE_ORDER == __LITTLE_ENDIAN
static inline uint64_t htonll(uint64_t x) { return bswap_64(x); }
static inline uint64_t ntohll(uint64_t x) { return bswap_64(x); }
#elif __BYTE_ORDER == __BIG_ENDIAN
static inline uint64_t htonll(uint64_t x) { return x; }
static inline uint64_t ntohll(uint64_t x) { return x; }
#else
#error __BYTE_ORDER is neither __LITTLE_ENDIAN nor __BIG_ENDIAN
#endif


#define NM (2+(1<<LM)) /* actual dimension including ghost cells for communications */
#define NV (ONE*(2+(1<<NDIM1))*(2+(1<<NDIM2))*(2+(1<<NDIM3))) /* size of rhs array */
#define NR (((NV+NM*NM+5*NM+7*LM+6)/7)*8) /* size of residual array */
#define MAXLEVEL (LT_DEFAULT+1) /* maximum number of levels */
#define M (NM+1) /* set at m=1024, can handle cases up to 1024^3 case */
#define MM (10)
#define	A (pow(5.0,13.0))
#define	X (314159265.0)
#define T_INIT 0
#define T_BENCH 1
#define T_MG3P 2
#define T_PSINV 3
#define T_RESID 4
#define T_RESID2 5
#define T_RPRJ3 6
#define T_INTERP 7
#define T_NORM2 8
#define T_COMM3 9
#define T_LAST 10



#if defined(USE_POW)
#define r23 pow(0.5, 23.0)
#define r46 (r23*r23)
#define t23 pow(2.0, 23.0)
#define t46 (t23*t23)
#else
#define r23 (0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5)
#define r46 (r23*r23)
#define t23 (2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0)
#define t46 (t23*t23)
#endif

/* global variables */
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
static int nx[MAXLEVEL+1];
static int ny[MAXLEVEL+1];
static int nz[MAXLEVEL+1];
static int m1[MAXLEVEL+1];
static int m2[MAXLEVEL+1];
static int m3[MAXLEVEL+1];
static int ir[MAXLEVEL+1];
static int debug_vec[8];
// static double u[NR];
// static double v[NV];
// static double r[NR];
#else
static int (*nx)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*ny)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*nz)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*m1)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*m2)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*m3)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*ir)=(int*)malloc(sizeof(int)*(MAXLEVEL+1));
static int (*debug_vec)=(int*)malloc(sizeof(int)*(8));

uint64_t u = 0;
uint64_t v = 0;
uint64_t r = 0;
//add a metadata .which range in the local, offset + length. check data on the remote or local.

//put them on the remote.
// static double (*u)=(double*)malloc(sizeof(double)*(NR));
// static double (*v)=(double*)malloc(sizeof(double)*(NV));
// static double (*r)=(double*)malloc(sizeof(double)*(NR));
#endif
static int is1, is2, is3, ie1, ie2, ie3, lt, lb;
static boolean timeron;

/* function prototypes */
static void bubble(double ten[][MM], int j1[][MM], int j2[][MM], int j3[][MM], int m, int ind);
static void comm3(struct resources *res, uint64_t u, int n1, int n2, int n3, int kk);
static void interp(struct resources *res, uint64_t z, int mm1, int mm2, int mm3, uint64_t u, int n1, int n2, int n3, int k);
static void mg3P(struct resources *res, uint64_t u, uint64_t v, uint64_t r, double a[4], double c[4], int n1, int n2, int n3, int k);
static void norm2u3(struct resources *res, uint64_t r, int n1, int n2, int n3, double* rnm2, double* rnmu, int nx, int ny, int nz);
static double power(double a, int n);
static void psinv(struct resources *res, uint64_t r, uint64_t u, int n1, int n2, int n3, double c[4], int k);
static void rep_nrm(struct resources *res, uint64_t u, int n1, int n2, int n3, char* title, int kk);
static void resid(struct resources *res, uint64_t u, uint64_t v, uint64_t r, int n1, int n2, int n3, double a[4], int k);
static void rprj3(struct resources *res, uint64_t r, int m1k, int m2k, int m3k, uint64_t s, int m1j, int m2j, int m3j, int k);
static void setup(int* n1, int* n2, int* n3, int k);
static void showall(struct resources *res, uint64_t z, int n1, int n2, int n3);
static void zero3(struct resources *res, uint64_t z, int n1, int n2, int n3);
static void zran3(struct resources *res, uint64_t z, int n1, int n2, int n3, int nx, int ny, int k);  
void vranlc(struct resources *res, int n, double *x_seed, double a, uint64_t y);
double randlc(double *x, double a);

/* mg */
int main(int argc, char *argv[]){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	/*
	 * -------------------------------------------------------------------------
	 * k is the current level. it is passed down through subroutine args
	 * and is not global. it is the current iteration
	 * -------------------------------------------------------------------------
	 */
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

	uint64_t addr = 0;
	//put variable u on the remote
	u_int64_t slices_cnt = ceil(sizeof(double) * NR / cpupool_mem_size) + 1;
	printf("slices is %d\n", slices_cnt);
	if (sizeof(double) * NR <= cpupool_mem_size ) {
		u = mymalloc(&res, sizeof(double) * NR, 0);
		if (poll_completion(&res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	} else {
		for (int i = 0; i < slices_cnt; i++) {
			if (i == slices_cnt - 1) {
				addr = mymalloc(&res, sizeof(double) * NR - (slices_cnt - 1) * cpupool_mem_size, i * cpupool_mem_size);
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				printf("memory limit is: %ld\n", addr + sizeof(double) * NR - (slices_cnt - 1) * cpupool_mem_size);
			} else {
				addr = mymalloc(&res, cpupool_mem_size, i * cpupool_mem_size); 
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				if (i == 0) {
					u = addr;
				}
			}
		}
	}
	printf("u: %ld u + NR %ld\n", u, u + NR * sizeof(double));

	// std::unordered_map<std::string, uint64_t> mymap;
	// mymap["u"] = u;
	// mymap["r"] = u + sizeof(double) * (NR);
	// mymap["v"] = u + sizeof(double) * (2 * NR);
	

	//put r on the remote
	addr = 0;
	slices_cnt = ceil(sizeof(double) * NR / cpupool_mem_size) + 1;
	printf("slices is %d\n", slices_cnt);
	if (sizeof(double) * NR <= cpupool_mem_size ) {
		r = mymalloc(&res, sizeof(double) * NR, sizeof(double) * (NR + 1000));
		if (poll_completion(&res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	} else {
		for (int i = 0; i < slices_cnt; i++) {
			if (i == slices_cnt - 1) {
				addr = mymalloc(&res, sizeof(double) * NR - (slices_cnt - 1) * cpupool_mem_size, sizeof(double) * (NR + 1000) + i * cpupool_mem_size);
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				printf("memory limit is: %ld\n", addr + sizeof(double) * NR - (slices_cnt - 1) * cpupool_mem_size);
			} else {
				addr = mymalloc(&res, cpupool_mem_size, sizeof(double) * (NR + 1000) + i * cpupool_mem_size); 
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				if (i == 0) {
					r = addr;
				}
			}
		}
	}
	printf("r: %ld r + NR %ld\n", r, r + NR * sizeof(double));
	//exit(0);

	addr = 0;
	//put variable v on the remote
	slices_cnt = ceil(sizeof(double) * NV / cpupool_mem_size) + 1;
	printf("slices is %d\n", slices_cnt);
	if (sizeof(double) * NV <= cpupool_mem_size ) {
		v = mymalloc(&res, sizeof(double) * NV, sizeof(double) * (2 * (NR + 1000)));
		if (poll_completion(&res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	} else {
		for (int i = 0; i < slices_cnt; i++) {
			if (i == slices_cnt - 1) {
				addr = mymalloc(&res, sizeof(double) * NV - (slices_cnt - 1) * cpupool_mem_size, sizeof(double) * (2 * (NR + 1000)) + i * cpupool_mem_size);
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				printf("memory limit is: %ld\n", addr + sizeof(double) * NV - (slices_cnt - 1) * cpupool_mem_size);
			} else {
				addr = mymalloc(&res, cpupool_mem_size, sizeof(double) * (2 * (NR + 1000)) + i * cpupool_mem_size); 
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				if (i == 0) {
					v = addr;
				}
			}
		}
	}
	printf("v: %ld v + NV %ld\n", v, v + NV * sizeof(double));

	

	int k, it;
	double t, tinit, mflops;

	double a[4], c[4];

	double rnm2, rnmu, epsilon;
	int n1, n2, n3, nit;
	double nn, verify_value, err;
	boolean verified;
	char class_npb;

	


	int i;
	char* t_names[T_LAST];
	double tmax;

	for(i=T_INIT; i<T_LAST; i++){
		timer_clear(i);
	}

	timer_start(T_INIT);	

	/*
	 * ----------------------------------------------------------------------
	 * read in and broadcast input data
	 * ----------------------------------------------------------------------
	 */
	FILE* fp;
	if((fp = fopen("timer.flag", "r")) != NULL){
		timeron = TRUE;
		t_names[T_INIT] = (char*) "init";
		t_names[T_BENCH] = (char*) "benchmk";
		t_names[T_MG3P] = (char*) "mg3P";
		t_names[T_PSINV] = (char*) "psinv";
		t_names[T_RESID] = (char*) "resid";
		t_names[T_RPRJ3] = (char*) "rprj3";
		t_names[T_INTERP] = (char*) "interp";
		t_names[T_NORM2] = (char*) "norm2";
		t_names[T_COMM3] = (char*) "comm3";
		fclose(fp);
	}else{
		timeron = FALSE;
	}
	fp = fopen("mg.input", "r");
	if(fp != NULL){
		printf(" Reading from input file mg.input\n");
		if(fscanf(fp, "%d", &lt) != 1){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		if(fscanf(fp, "%d%d%d", &nx[lt], &ny[lt], &nz[lt]) != 3){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		if(fscanf(fp, "%d", &nit) != 1){
			printf(" Error in reading elements\n");
			exit(1);
		}
		while(fgetc(fp) != '\n');
		for(i = 0; i <= 7; i++) {
			if(fscanf(fp, "%d", &debug_vec[i]) != 1){
				printf(" Error in reading elements\n");
				exit(1);
			}
		}
		fclose(fp);
	}else{
		printf(" No input file. Using compiled defaults\n");
		lt = LT_DEFAULT;
		nit = NIT_DEFAULT;
		nx[lt] = NX_DEFAULT;
		ny[lt] = NY_DEFAULT;
		nz[lt] = NZ_DEFAULT;
		for(i = 0; i <= 7; i++){
			debug_vec[i] = DEBUG_DEFAULT;
		}
		//debug_vec[2] = 2;
	}

	if((nx[lt] != ny[lt]) || (nx[lt] != nz[lt])){
		class_npb = 'U';
	}else if(nx[lt] == 32 && nit == 4){
		class_npb = 'S';
	}else if(nx[lt] == 128 && nit == 4){
		class_npb = 'W';
	}else if(nx[lt] == 256 && nit == 4){
		class_npb = 'A';
	}else if(nx[lt] == 256 && nit == 20){
		class_npb = 'B';
	}else if(nx[lt] == 512 && nit == 20){
		class_npb = 'C';
	}else if(nx[lt] == 1024 && nit == 50){  
		class_npb = 'D';
	}else if(nx[lt] == 2048 && nit == 50){  
		class_npb = 'E';	
	}else{
		class_npb = 'U';
	}

	/*
	 * ---------------------------------------------------------------------
	 * use these for debug info:
	 * ---------------------------------------------------------------------
	 * debug_vec(0) = 1 !=> report all norms
	 * debug_vec(1) = 1 !=> some setup information
	 * debug_vec(1) = 2 !=> more setup information
	 * debug_vec(2) = k => at level k or below, show result of resid
	 * debug_vec(3) = k => at level k or below, show result of psinv
	 * debug_vec(4) = k => at level k or below, show result of rprj
	 * debug_vec(5) = k => at level k or below, show result of interp
	 * debug_vec(6) = 1 => (unused)
	 * debug_vec(7) = 1 => (unused)
	 * ---------------------------------------------------------------------
	 */
	a[0] = -8.0/3.0;
	a[1] =  0.0;
	a[2] =  1.0/6.0;
	a[3] =  1.0/12.0;

	if(class_npb == 'A' || class_npb == 'S' || class_npb =='W'){
		/* coefficients for the s(a) smoother */
		c[0] =  -3.0/8.0;
		c[1] =  +1.0/32.0;
		c[2] =  -1.0/64.0;
		c[3] =   0.0;
	}else{
		/* coefficients for the s(b) smoother */
		c[0] =  -3.0/17.0;
		c[1] =  +1.0/33.0;
		c[2] =  -1.0/61.0;
		c[3] =   0.0;
	}

	lb = 1;
	k = lt;

	setup(&n1,&n2,&n3,k);
	
	zero3(&res, u,n1,n2,n3);
	
	// int cnt = ceil(NR * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, u + i * cpupool_mem_size, res.buf, NR * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, u + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int i = 0; i < min(NR * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); i++) {
	// 		printf("%f \n", res.buf[i]);
	// 	}
	// }
	// exit(0); 

	zran3(&res, v, n1,n2,n3,nx[lt],ny[lt],k);	
	

	
	norm2u3(&res, v,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
	// int cnt = ceil(NV * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, v + i * cpupool_mem_size, res.buf, NV * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, v + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int i = 0; i < min(NV * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); i++) {
	// 		printf("%f \n", res.buf[i]);
	// 	}
	// }
	// exit(0); 
	// printf("rnm2 = %f, rnmu = %f", rnm2, rnmu);
	// exit(0);
	

	printf("\n\n NAS Parallel Benchmarks 4.1 Parallel C++ version with OpenMP - MG Benchmark\n\n");
	printf(" Size: %3dx%3dx%3d (class_npb %1c)\n", nx[lt], ny[lt], nz[lt], class_npb);
	printf(" Iterations: %3d\n", nit);
	

	#pragma omp parallel
	{
		resid(&res,u,v,r,n1,n2,n3,a,k);
		norm2u3(&res, r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);	
		mg3P(&res, u,v,r,a,c,n1,n2,n3,k);
		resid(&res, u,v,r,n1,n2,n3,a,k);
	}
	
	setup(&n1,&n2,&n3,k);

	zero3(&res, u,n1,n2,n3);

	zran3(&res, v,n1,n2,n3,nx[lt],ny[lt],k);
	timer_stop(T_INIT);
	tinit = timer_read(T_INIT);
	printf(" Initialization time: %15.3f seconds\n", tinit);

	for(i=T_BENCH; i<T_LAST; i++){
		timer_clear(i);
	}
	timer_start(T_BENCH);

	#pragma omp parallel firstprivate(nit) private(it)
    {

		if(timeron){
			#pragma omp master
				timer_start(T_RESID2);
		}
		
		resid(&res, u,v,r,n1,n2,n3,a,k);
		if(timeron){
			#pragma omp master
				timer_stop(T_RESID2);
		}

		norm2u3(&res, r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
		for(it = 1; it <= nit; it++){
			if((it==1)||(it==nit)||((it%5)==0)){
				#pragma omp master
					printf("  iter %3d\n",it);
			}
			
			if(timeron){
				#pragma omp master
					timer_start(T_MG3P);
			}
			mg3P(&res, u,v,r,a,c,n1,n2,n3,k);
			if(timeron){
				#pragma omp master
					timer_stop(T_MG3P);
			}
			if(timeron){
				#pragma omp master
					timer_start(T_RESID2);
			}

			resid(&res, u,v,r,n1,n2,n3,a,k);
	
			if(timeron){
				#pragma omp master
					timer_stop(T_RESID2);
			}
		}
	
		norm2u3(&res,r,n1,n2,n3,&rnm2,&rnmu,nx[lt],ny[lt],nz[lt]);
		
	} /* end parallel */

	timer_stop(T_BENCH);
	// FILE *file_v, *file_r, *file_u;
    // file_v = fopen("v.txt", "w"); 
	// file_r = fopen("r.txt", "w");
	// file_u = fopen("u.txt", "w");

    // if (file_v == NULL) {
    //     perror("Error opening file");
    // }
	// if (file_u == NULL) {
    //     perror("Error opening file");
    // }
	// if (file_r == NULL) {
    //     perror("Error opening file");
    // }
	
	// int cnt = ceil(NR * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, u + i * cpupool_mem_size, res.buf, NR * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, u + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NR * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_u, "%f \n", res.buf[j]);
	// 	}
	// }
	
	
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, r + i * cpupool_mem_size, res.buf, NR * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, r + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NR * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_r, "%f \n", res.buf[j]);
	// 	}
	// }

	// cnt = ceil(NV * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, v + i * cpupool_mem_size, res.buf, NV * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, v + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NV * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_v, "%f \n", res.buf[j]);
	// 	}
	// }
	// fclose(file_u);
	// fclose(file_v);
	// fclose(file_r);
	// exit(0);
	t = timer_read(T_BENCH);    	

	verified = FALSE;
	verify_value = 0.0;	

	printf(" Benchmark completed\n");

	epsilon = 1.0e-8;
	if(class_npb != 'U'){
		if(class_npb == 'S'){
			verify_value = 0.5307707005734e-04;
		}else if(class_npb == 'W'){
			verify_value = 0.6467329375339e-05;
		}else if(class_npb == 'A'){
			verify_value = 0.2433365309069e-05;
		}else if(class_npb == 'B'){
			verify_value = 0.1800564401355e-05;
		}else if(class_npb == 'C'){
			verify_value = 0.5706732285740e-06;
		}else if(class_npb == 'D'){
			verify_value = 0.1583275060440e-09;
		}else if(class_npb == 'E'){
			verify_value = 0.8157592357404e-10; 
		}

		err = fabs(rnm2-verify_value) / verify_value;
		if(err <= epsilon){
			verified = TRUE;
			printf(" VERIFICATION SUCCESSFUL\n");
			printf(" L2 Norm is %20.13e\n", rnm2);
			printf(" Error is   %20.13e\n", err);
		}else{
			verified = FALSE;
			printf(" VERIFICATION FAILED\n");
			printf(" L2 Norm is             %20.13e\n", rnm2);
			printf(" The correct L2 Norm is %20.13e\n", verify_value);
		}
	}else{
		verified = FALSE;
		printf(" Problem size unknown\n");
		printf(" NO VERIFICATION PERFORMED\n");
	}

	nn = 1.0*nx[lt]*ny[lt]*nz[lt];

	if(t!=0.0){
		mflops = 58.0*nit*nn*1.0e-6/t;
	}else{
		mflops = 0.0;
	}

	setenv("OMP_NUM_THREADS","1",0);
	c_print_results((char*)"MG",
			class_npb,
			nx[lt],
			ny[lt],
			nz[lt],
			nit,
			t,
			mflops,
			(char*)"          floating point",
			verified,
			(char*)NPBVERSION,
			(char*)COMPILETIME,
			(char*)COMPILERVERSION,
			(char*)LIBVERSION,
			std::getenv("OMP_NUM_THREADS"),
			(char*)CS1,
			(char*)CS2,
			(char*)CS3,
			(char*)CS4,
			(char*)CS5,
			(char*)CS6,
			(char*)CS7);

	/*
	 * ---------------------------------------------------------------------
	 * more timers
	 * ---------------------------------------------------------------------
	 */
	if(timeron){
		tmax = timer_read(T_BENCH);
		if(tmax==0.0){tmax=1.0;}
		printf("  SECTION   Time (secs)\n");
		for(i=T_BENCH; i<T_LAST; i++){
			t = timer_read(i);
			if(i==T_RESID2){
				t = timer_read(T_RESID) - t;
				printf("    --> %8s:%9.3f  (%6.2f%%)\n", "mg-resid", t, t*100.0/tmax);
			}else{
				printf("  %-8s:%9.3f  (%6.2f%%)\n", t_names[i], t, t*100.0/tmax);
			}
		}
	}

	return 0;
}

/*
 * ---------------------------------------------------------------------
 * bubble does a bubble sort in direction dir
 * ---------------------------------------------------------------------
 */
static void bubble(double ten[][MM], int j1[][MM], int j2[][MM], int j3[][MM], int m, int ind){
	double temp;
	int i, j_temp;

	if(ind == 1){
		for(i = 0; i < m-1; i++){
			if(ten[ind][i] > ten[ind][i+1]){
				temp = ten[ind][i+1];
				ten[ind][i+1] = ten[ind][i];
				ten[ind][i] = temp;

				j_temp = j1[ind][i+1];
				j1[ind][i+1] = j1[ind][i];
				j1[ind][i] = j_temp;

				j_temp = j2[ind][i+1];
				j2[ind][i+1] = j2[ind][i];
				j2[ind][i] = j_temp;

				j_temp = j3[ind][i+1];
				j3[ind][i+1] = j3[ind][i];
				j3[ind][i] = j_temp;
			}else{
				return;
			}
		}
	}else{
		for(i = 0; i < m-1; i++){
			if(ten[ind][i] < ten[ind][i+1]){
				temp = ten[ind][i+1];
				ten[ind][i+1] = ten[ind][i];
				ten[ind][i] = temp;

				j_temp = j1[ind][i+1];
				j1[ind][i+1] = j1[ind][i];
				j1[ind][i] = j_temp;

				j_temp = j2[ind][i+1];
				j2[ind][i+1] = j2[ind][i];
				j2[ind][i] = j_temp;

				j_temp = j3[ind][i+1];
				j3[ind][i+1] = j3[ind][i];
				j3[ind][i] = j_temp;
			}else{
				return;
			}
		}
	}
}

/*
 * ---------------------------------------------------------------------
 * comm3 organizes the communication on all borders 
 * ---------------------------------------------------------------------
 */
static void comm3(struct resources *res, uint64_t u, int n1, int n2, int n3, int kk){
// #ifdef __clang__
// 		using custom_cast = double (*)[n2][n1];
// 		custom_cast u = reinterpret_cast<custom_cast>(pointer_u);
// #else
// 		double (*u)[n2][n1] = (double (*)[n2][n1])pointer_u;
// #endif
	// struct timespec req, rem;
	// req.tv_sec = 0;       
    // req.tv_nsec = 100; 
	int i1, i2, i3;
	if(timeron){
		#pragma omp master
			timer_start(T_COMM3);
	}
	#pragma omp for
		/* axis = 1 */
	for(i3 = 1; i3 < n3-1; i3++){
		//int tid = omp_get_thread_num();
		//int num_element = min(n2 * n1, static_cast<int>(block_size / sizeof(double)));	
		//for (int i = 0; i < ceil(n2 * n1 / num_element) + 1; i++) {
			myread(res, u + (i3 * n2 * n1) * sizeof(double), res->buf, sizeof(double) * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}

			for(i2 = 1; i2 < n2 - 1; i2++){
				mywrite(res, u + (i3 * n2 * n1 + i2 * n1) * sizeof(double), res->buf + i2 * n1 + n1 - 2, sizeof(double));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
	
				mywrite(res, u + (i3 * n2 * n1 + i2 * n1 + n1 - 1) * sizeof(double), res->buf + i2 * n1 + 1, sizeof(double));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
			//******************************
			//u[i3][i2][0] = u[i3][i2][n1-2];
			//u[i3][i2][n1-1] = u[i3][i2][1];			
			}
		//}
		/* axis = 2 */
		//for (int i = 0; i < ceil(n2 * n1 / num_element) + 1; i++) {
			// myread(res, u + (i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * n2 * n1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }

			for (int i1 = 0; i1 < n1; i1++) {
				mywrite(res, u + (i3 * n2 * n1 + i1) * sizeof(double), res->buf + (n2 - 2) * n1 + i1, sizeof(double));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
	
				mywrite(res, u + (i3 * n2 * n1 + (n2 - 1) * n1 + i1) * sizeof(double), res->buf + 1 * n1 + i1, sizeof(double));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
			}
			
			// u[i3][0][i1] = u[i3][n2-2][i1];
			// u[i3][n2-1][i1] = u[i3][1][i1];			
		//}
	}

	


	//int num_element = min(n2 * n1, static_cast<int>(block_size / sizeof(double)));
	/* axis = 3 */

 	//u[0][i2][i1] = u[n3-2][i2][i1];
	// 	//u[n3-1][i2][i1] = u[1][i2][i1];	
	//#pragma omp for
	//for (int i = 0; i < 1 ; i++) {
		//int tid = omp_get_thread_num();
		myread(res, u + (n3 - 2) * n2 * n1 * sizeof(double), res->buf, sizeof(double) * n2 * n1);
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}

		myread(res, u + (1 * n2 * n1 ) * sizeof(double), res->buf + n2 * n1, sizeof(double) * n2 * n1);
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}

		mywrite(res, u, res->buf, sizeof(double) * n2 * n1);
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
		
		mywrite(res, u + ((n3 - 1) * n2 * n1) * sizeof(double) , res->buf + n2 * n1, sizeof(double) * n2 * n1);
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	//}
	// #pragma omp for
	// for(i2 = 0; i2 < n2; i2++){
	// 	for(i1 = 0; i1 < n1; i1++){
	// 		myread(res, u + ((n3 - 2) * n2 * n1 + i2 * n1 + i1) * sizeof(double), res->buf, sizeof(double));
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 		mywrite(res, u + (i2 * n1 + i1) * sizeof(double), res->buf, sizeof(double));
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 		myread(res, u + (n2 * n1 + i2 * n1 + i1) * sizeof(double), res->buf, sizeof(double));
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 		mywrite(res, u + ((n3 - 1) * n2 * n1 + i2 * n1 + i1) * sizeof(double), res->buf, sizeof(double));
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
			
	// 	}
	// }
		// myread(res, u + (n2 * n1) * sizeof(double), res->buf + block_size / sizeof(double), sizeof(double) * n2 * n1);
		// if (poll_completion(res)) {
		// 	fprintf(stderr, "poll completion failed\n");
		// }	
		// mywrite(res, u + ((n3 - 1) * n2 * n1) * sizeof(double), res->buf + block_size / sizeof(double), sizeof(double) * n2 * n1);
		// if (poll_completion(res)) {
		// 	fprintf(stderr, "poll completion failed\n");
		// }
		
	
	// for(i2 = 0; i2 < n2; i2++){
	// 	for(i1 = 0; i1 < n1; i1++){
	// 	//u[0][i2][i1] = u[n3-2][i2][i1];
	// 	//u[n3-1][i2][i1] = u[1][i2][i1];			
	// 	}
	// }

	if(timeron){
		#pragma omp master
			timer_stop(T_COMM3);
	}
}


/*
 * --------------------------------------------------------------------
 * interp adds the trilinear interpolation of the correction
 * from the coarser grid to the current approximation: u = u + Qu'
 *     
 * observe that this  implementation costs  16A + 4M, where
 * A and M denote the costs of addition and multiplication.  
 * note that this vectorizes, and is also fine for cache 
 * based machines. vector machines may get slightly better 
 * performance however, with 8 separate "do i1" loops, rather than 4.
 * --------------------------------------------------------------------
 */
static void 
interp(struct resources *res, uint64_t z, int mm1, int mm2, int mm3, uint64_t u, int n1, int n2, int n3, int k){
// #ifdef __clang__
// 	using custom_cast = double (*)[mm2][mm1];
// 	custom_cast z = reinterpret_cast<custom_cast>(pointer_z);
// 	using custom_cast2 = double (*)[n2][n1];
// 	custom_cast2 u = reinterpret_cast<custom_cast2>(pointer_u);
// #else
// 	double (*z)[mm2][mm1] = (double (*)[mm2][mm1])pointer_z;
// 	double (*u)[n2][n1] = (double (*)[n2][n1])pointer_u;
// #endif

	int i3, i2, i1, d1, d2, d3, t1, t2, t3;

	/* 
	 * --------------------------------------------------------------------
	 * note that m = 1037 in globals.h but for this only need to be
	 * 535 to handle up to 1024^3
	 * integer m
	 * parameter( m=535 )
	 * --------------------------------------------------------------------
	 */
	double z1[M], z2[M], z3[M];

	if(timeron){
		#pragma omp master
			timer_start(T_INTERP);
	}

	if(n1 != 3 && n2 != 3 && n3 != 3){
		#pragma omp for
		for(i3 = 0; i3 < mm3-1; i3++){
			int tid = 0; 
			int num_element = min(2 * mm2 * mm1 + 2 * n1 * n2, static_cast<int>(block_size / sizeof(double)));
				for (int i = 0; i < 1; i++) {
					myread(res, z + (i3 * mm2 * mm1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					myread(res, u + (2 * i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				for(i2 = 0; i2 < mm2-1; i2++){

				for(i1 = 0; i1 < mm1; i1++){
					z1[i1] = res->buf[tid * block_size / sizeof(double) + (i2 + 1) * mm1 + i1] + res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1];
					z2[i1] = res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + i2 * mm1 + i1] + res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1];
					z3[i1] = res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 + 1) * mm1 + i1] + res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + i2 * mm1 + i1] + z1[i1];
					// z1[i1] = z[i3][i2+1][i1] + z[i3][i2][i1];
					// z2[i1] = z[i3+1][i2][i1] + z[i3][i2][i1];
					// z3[i1] = z[i3+1][i2+1][i1] + z[i3+1][i2][i1] + z1[i1];
					//printf("%f %f %f \n",z1[i1], z2[i1], z3[i1]);
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1] += res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1];
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1 + 1] +=  0.5 * (res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1 + 1] + res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1]);
					mywrite(res, u + (2 * i3 * n2 * n1 + 2 * i2 * n1 +2 * i1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1, sizeof(double) * 2);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					// u[2*i3][2*i2][2*i1] = u[2*i3][2*i2][2*i1] + z[i3][i2][i1];
					// u[2*i3][2*i2][2*i1+1] = u[2*i3][2*i2][2*i1+1] + 0.5 * (z[i3][i2][i1+1]+z[i3][i2][i1]);
				}
				
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1] += 0.5 * z1[i1];
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.25 * (z1[i1] + z1[i1 + 1]);
					mywrite(res, u + (2 * i3 * n2 * n1 + (2 * i2 + 1) * n1 + 2 * i1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1, sizeof(double) * 2);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					// u[2*i3][2*i2+1][2*i1] = u[2*i3][2*i2+1][2*i1] + 0.5 * z1[i1];
					// u[2*i3][2*i2+1][2*i1+1] = u[2*i3][2*i2+1][2*i1+1] + 0.25 * ( z1[i1] + z1[i1+1] );
					
				}

				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1] += 0.5 * z2[i1];
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1 + 1] += 0.25 * (z2[i1] + z2[i1 + 1]);
					mywrite(res, u + ((2 * i3 + 1) * n2 * n1 + (2 * i2) * n1 + 2 * i1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1, sizeof(double) * 2);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					// u[2*i3+1][2*i2][2*i1] = u[2*i3+1][2*i2][2*i1] + 0.5 * z2[i1];
					// u[2*i3+1][2*i2][2*i1+1] = u[2*i3+1][2*i2][2*i1+1]+ 0.25*( z2[i1] + z2[i1+1] );
					
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1] += 0.25 * z3[i1];
					res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.125 * (z3[i1] + z3[i1 + 1]);
					mywrite(res, u + ((2 * i3 + 1) * n2 * n1 + (2 * i2 + 1) * n1 + 2 * i1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1, sizeof(double) * 2);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					// u[2*i3+1][2*i2+1][2*i1] = u[2*i3+1][2*i2+1][2*i1] + 0.25* z3[i1];
					// u[2*i3+1][2*i2+1][2*i1+1] = u[2*i3+1][2*i2+1][2*i1+1] + 0.125*( z3[i1] + z3[i1+1] );
					//printf("%f %f\n", res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1], res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1 + 1]);
				}
			}	
			
			// mywrite(res, z + (i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
			// mywrite(res, u + (2 * i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
		}
	 
			// for(i2 = 0; i2 < mm2-1; i2++){

			// 	for(i1 = 0; i1 < mm1; i1++){
			// 		z1[i1] = res.buf[]
			// 		// z1[i1] = z[i3][i2+1][i1] + z[i3][i2][i1];
			// 		// z2[i1] = z[i3+1][i2][i1] + z[i3][i2][i1];
			// 		// z3[i1] = z[i3+1][i2+1][i1] + z[i3+1][i2][i1] + z1[i1];
			// 	}
			// 	for(i1 = 0; i1 < mm1-1; i1++){
			// 		u[2*i3][2*i2][2*i1] = u[2*i3][2*i2][2*i1]
			// 			+z[i3][i2][i1];
			// 		u[2*i3][2*i2][2*i1+1] = u[2*i3][2*i2][2*i1+1]
			// 			+0.5*(z[i3][i2][i1+1]+z[i3][i2][i1]);
			// 	}
			// 	for(i1 = 0; i1 < mm1-1; i1++){
			// 		u[2*i3][2*i2+1][2*i1] = u[2*i3][2*i2+1][2*i1]
			// 			+0.5 * z1[i1];
			// 		u[2*i3][2*i2+1][2*i1+1] = u[2*i3][2*i2+1][2*i1+1]
			// 			+0.25*( z1[i1] + z1[i1+1] );
			// 	}
			// 	for(i1 = 0; i1 < mm1-1; i1++){
			// 		u[2*i3+1][2*i2][2*i1] = u[2*i3+1][2*i2][2*i1]
			// 			+0.5 * z2[i1];
			// 		u[2*i3+1][2*i2][2*i1+1] = u[2*i3+1][2*i2][2*i1+1]
			// 			+0.25*( z2[i1] + z2[i1+1] );
			// 	}
			// 	for(i1 = 0; i1 < mm1-1; i1++){
			// 		u[2*i3+1][2*i2+1][2*i1] = u[2*i3+1][2*i2+1][2*i1]
			// 			+0.25* z3[i1];
			// 		u[2*i3+1][2*i2+1][2*i1+1] = u[2*i3+1][2*i2+1][2*i1+1]
			// 			+0.125*( z3[i1] + z3[i1+1] );
			// 	}
			// }
		}

	}else{
		if(n1 == 3){
			d1 = 2;
			t1 = 1;
		}else{
			d1 = 1;
			t1 = 0;
		}      
		if(n2 == 3){
			d2 = 2;
			t2 = 1;
		}else{
			d2 = 1;
			t2 = 0;
		}          
		if(n3 == 3){
			d3 = 2;
			t3 = 1;
		}else{
			d3 = 1;
			t3 = 0;
		}
		#pragma omp for
		for(i3 = d3; i3 <= mm3-1; i3++){
			int tid = 0; 
			int num_element = min(2 * n2 * n1, static_cast<int>(block_size / sizeof(double)));
				for (int i = 0; i < 0 + 1; i++) {
					myread(res, z + ((i3 - 1) * mm2 * mm1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * mm2 * mm1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}

					myread(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1, sizeof(double) *  n2 * n1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					for(i2 = d2; i2 <= mm2-1; i2++){
						for(i1 = d1; i1 <= mm1-1; i1++){
							res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - d2 -1) * n1 + 2 * i1 - d1 - 1] += res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1 -1];
							mywrite(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - d1 - 1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - d2 -1) * n1 + 2 * i1 - d1 - 1, sizeof(double));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
							//u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] = u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] + z[i3-1][i2-1][i1-1];
						}
						for(i1 = 1; i1 <= mm1-1; i1++){
							res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - d2 -1) * n1 + 2 * i1 - t1 - 1] += 0.5 * (res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1] 
																																		+ res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1 -1] );
							mywrite(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - t1 - 1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - d2 -1) * n1 + 2 * i1 - t1 - 1, sizeof(double));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
							//u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] = u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] + 0.5*(z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
						}
					}
					for(i2 = 1; i2 <= mm2-1; i2++){
						for (i1 = d1; i1 <= mm1-1; i1++) {
							res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - d1 - 1] += 0.5 * (res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1 -1 ] 
																																+ res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1 - 1]);
							mywrite(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - d1 - 1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - d1 - 1, sizeof(double));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
							//u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] = u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] + 0.5*(z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
						}
						for(i1 = 1; i1 <= mm1-1; i1++){
							res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - t1 - 1] += 0.25 * (res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1] +
																																res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1] +
																																res->buf[tid * block_size / sizeof(double) + i2 * mm1 + i1 - 1] + 
																																res->buf[tid * block_size / sizeof(double) + (i2 - 1) * mm1 + i1 -1]);
							mywrite(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - t1 - 1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1 + (2 * i2 - t2 - 1) * n1 + 2 * i1 - t1 - 1, sizeof(double));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
							//u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] = u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] + 0.25*(z[i3-1][i2][i1]+z[i3-1][i2-1][i1] + z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
						}
					}
					// mywrite(res, z + ((i3 - 1) * mm2 * mm1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * mm2 * mm1);
					// if (poll_completion(res)) {
					// 	fprintf(stderr, "poll completion failed\n");
					// }
					// mywrite(res, u + ((2 * i3 - d3 - 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + mm2 * mm1, sizeof(double) * n2 * n1);
					// if (poll_completion(res)) {
					// 	fprintf(stderr, "poll completion failed\n");
					// }
				}
			
			// for(i2 = d2; i2 <= mm2-1; i2++){
			// 	for(i1 = d1; i1 <= mm1-1; i1++){
			// 		u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] = u[2*i3-d3-1][2*i2-d2-1][2*i1-d1-1] + z[i3-1][i2-1][i1-1];
			// 	}
			// 	for(i1 = 1; i1 <= mm1-1; i1++){
			// 		u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] = u[2*i3-d3-1][2*i2-d2-1][2*i1-t1-1] + 0.5*(z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
			// 	}
			// }
			// for(i2 = 1; i2 <= mm2-1; i2++){
			// 	for ( i1 = d1; i1 <= mm1-1; i1++) {
			// 		u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] = u[2*i3-d3-1][2*i2-t2-1][2*i1-d1-1] + 0.5*(z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
			// 	}
			// 	for(i1 = 1; i1 <= mm1-1; i1++){
			// 		u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] = u[2*i3-d3-1][2*i2-t2-1][2*i1-t1-1] + 0.25*(z[i3-1][i2][i1]+z[i3-1][i2-1][i1] + z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
			// 	}
			// }
		}
		#pragma omp for
		for(i3 = 1; i3 <= mm3-1; i3++){
			int tid = 0; 
			int num_element = min(n2 * n1, static_cast<int>(block_size / sizeof(double)));
				for (int i = 0; i < 0 + 1; i++) {
					myread(res, z + ((i3 - 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
					myread(res, u + ((2 * i3 - t3 - 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * n2 * n1);
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
			for(i2 = d2; i2 <= mm2-1; i2++){
				for(i1 = d1; i1 <= mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) +  2 * mm2 * mm1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - d1 - 1] += 0.5 * (res->buf[tid * block_size / sizeof(double) + n2 * n1 + (i2 - 1) * n1 + i1 -1] 
																													+ res->buf[tid * block_size  / sizeof(double) + (i2 - 1) * n1 + i1 -1]);
					//u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1] = u[2*i3-t3-1][2*i2-d2-1][2*i1-d1-1] +0.5*(z[i3][i2-1][i1-1]+z[i3-1][i2-1][i1-1]);
					mywrite(res, u + ((2 * i3 - t3 - 1) *  n2 * n1 + (2 * i2 -d2 -1) * n1 + 2*i1-d1-1) * sizeof(double), res->buf + tid * block_size / sizeof(double) +  2 * mm2 * mm1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - d1 - 1, sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}
				for(i1 = 1; i1 <= mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) +  2 * mm2 * mm1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - t1 - 1] += 0.25 * (res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 - 1) * n1 + i1] + 
																															res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 - 1) * n1 + i1 -1]+
																															res->buf[tid * block_size  / sizeof(double)+ (i2 - 1) * n1 + i1] + 
																															res->buf[tid * block_size  / sizeof(double)+ (i2 - 1) * n1 + i1 -1]);
					// u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1] = u[2*i3-t3-1][2*i2-d2-1][2*i1-t1-1]
					// +0.25*(z[i3][i2-1][i1]+z[i3][i2-1][i1-1]
					// +z[i3-1][i2-1][i1]+z[i3-1][i2-1][i1-1]);
					mywrite(res, u + ((2 * i3 - t3 - 1) *  n2 * n1 + (2 * i2 -d2 -1) * n1 + 2*i1-t1-1) * sizeof(double), res->buf + tid * block_size / sizeof(double) +  2 * mm2 * mm1 + (2 * i2 - d2 - 1) * n1 + 2 * i1 - t1 - 1, sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}
			}
			for(i2 = 1; i2 <= mm2-1; i2++){
				for (i1 = d1; i1 <= mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) +  2 * n2 * n1 + (2*i2-t2-1) * n1 + 2*i1-d1-1] += 0.25 * (res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + i2 * n1 + i1 - 1] + 
																															res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 - 1) * n1 + i1 -1]+
																															res->buf[tid * block_size  / sizeof(double)+ i2 * n1 + i1 -1 ] + 
																															res->buf[tid * block_size  / sizeof(double)+ (i2 - 1) * n1 + i1 -1]);
					// u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1] = u[2*i3-t3-1][2*i2-t2-1][2*i1-d1-1]+
					// 0.25*(z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
					// +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
					mywrite(res, u + ((2 * i3 - t3 - 1) *  n2 * n1 + (2 * i2 -t2 -1) * n1 + 2*i1-d1-1) * sizeof(double), res->buf + tid * block_size / sizeof(double) +  2 * n2 * n1 + (2*i2-t2-1) * n1 + 2*i1-d1-1, sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}
				for(i1 = 1; i1 <= mm1-1; i1++){
					res->buf[tid * block_size / sizeof(double) +  2 * n2 * n1 + (2*i2-t2-1) * n1 + 2*i1-t1-1] += 0.125 * (res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + i2 * n1 + i1] + 
																															res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 - 1) * n1 + i1]+
																															res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + i2 * n1 + i1 - 1]+
																															res->buf[tid * block_size / sizeof(double) + mm2 * mm1 + (i2 - 1) * n1 + i1 - 1]+
																															res->buf[tid * block_size  / sizeof(double)+ i2 * n1 + i1] + 
																															res->buf[tid * block_size  / sizeof(double)+ (i2 - 1) * n1 + i1] +
																															res->buf[tid * block_size  / sizeof(double)+ i2 * n1 + i1 -1] + 
																															res->buf[tid * block_size  / sizeof(double)+ (i2 - 1) * n1 + i1 - 1]);
																															
					// u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1] = u[2*i3-t3-1][2*i2-t2-1][2*i1-t1-1] 
					// +0.125*(z[i3][i2][i1]+z[i3][i2-1][i1]
					// +z[i3][i2][i1-1]+z[i3][i2-1][i1-1]
					// +z[i3-1][i2][i1]+z[i3-1][i2-1][i1]
					// +z[i3-1][i2][i1-1]+z[i3-1][i2-1][i1-1]);
					mywrite(res, u + ((2 * i3 - t3 - 1) *  n2 * n1 + (2 * i2 -t2 -1) * n1 + 2*i1-t1-1) * sizeof(double), res->buf + tid * block_size / sizeof(double) +  2 * n2 * n1 + (2*i2-t2-1) * n1 + 2*i1-t1-1, sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}
			}
			// mywrite(res, z + ((i3 - 1) * n2 * n1 + i * n2 * n1) * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
			// mywrite(res, u + ((2 * i3 - t3 - 1) *  n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * n2 * n1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
		}
	}
	if(timeron){
		#pragma omp master
			timer_stop(T_INTERP);
	}
	#pragma omp single
    {
		if(debug_vec[0] >= 1){
			rep_nrm(res,z,mm1,mm2,mm3,(char*)"z: inter",k-1);
			rep_nrm(res,u,n1,n2,n3,(char*)"u: inter",k);
		}
		if(debug_vec[5] >= k){
			showall(res,z,mm1,mm2,mm3);
			showall(res,u,n1,n2,n3);
		}
	}
}
}

/* 
 * --------------------------------------------------------------------
 * multigrid v-cycle routine
 * --------------------------------------------------------------------
 */
static void mg3P(struct resources *res, uint64_t u, uint64_t v, uint64_t r, double a[4], double c[4], int n1, int n2, int n3, int k){
	int j;

	/*
	 * --------------------------------------------------------------------
	 * down cycle.
	 * restrict the residual from the find grid to the coarse
	 * -------------------------------------------------------------------
	 */
	// myread(res, r, res->buf, sizeof(double));
	// if (poll_completion(res)) {
	// 	fprintf(stderr, "poll completion failed\n");
	// }
	// myread(res, r + 39304 * sizeof(double), res->buf + sizeof(double) * 1, sizeof(double));
	// if (poll_completion(res)) {
	// 	fprintf(stderr, "poll completion failed\n");
	// }
	// myread(res, r + 45136 * sizeof(double), res->buf + sizeof(double) * 2, sizeof(double));
	// if (poll_completion(res)) {
	// 	fprintf(stderr, "poll completion failed\n");
	// }
	// myread(res, r + 46136 * sizeof(double), res->buf + sizeof(double) * 3, sizeof(double));
	// if (poll_completion(res)) {
	// 	fprintf(stderr, "poll completion failed\n");
	// }
	// myread(res, r + 46352 * sizeof(double), res->buf + sizeof(double) * 4, sizeof(double));
	// if (poll_completion(res)) {
	// 	fprintf(stderr, "poll completion failed\n");
	// }
	// printf("%f, %f, %f, %f, %f\n", res->buf[0], res->buf[1], res->buf[2],res->buf[3], res->buf[4]);
	for(k = lt; k >= lb+1; k--){
		j = k-1;
		rprj3(res, r + ir[k] * sizeof(double), m1[k], m2[k], m3[k], r + ir[j] * sizeof(double), m1[j], m2[j], m3[j], k);
		// myread(res, r + ir[k] * sizeof(double), res->buf, sizeof(double));
		// if (poll_completion(res)) {
		// 	fprintf(stderr, "poll completion failed\n");
		// }
		// myread(res, r + ir[j] * sizeof(double), res->buf + sizeof(double), sizeof(double));
		// if (poll_completion(res)) {
		// 	fprintf(stderr, "poll completion failed\n");
		// }
		// printf("%f, %f\n", res->buf[0], res->buf[1]);
		// printf("%d, %d\n", ir[k], ir[j]);
		// printf("%d %d %d %d %d %d \n\n", m1[k], m2[k], m3[k], m1[j], m2[j], m3[j]);
		//rprj3(res, &r[ir[k]], m1[k], m2[k], m3[k], &r[ir[j]], m1[j], m2[j], m3[j], k);
	}
	k = lb;
	/*
	 * --------------------------------------------------------------------
	 * compute an approximate solution on the coarsest grid
	 * --------------------------------------------------------------------
	 */
	zero3(res, u + ir[k] * sizeof(double), m1[k], m2[k], m3[k]);
	
	psinv(res, r + ir[k] * sizeof(double), u + ir[k] * sizeof(double), m1[k], m2[k], m3[k], c, k);
	// zero3(res, &u[ir[k]], m1[k], m2[k], m3[k]);
	// psinv(res, &r[ir[k]], &u[ir[k]], m1[k], m2[k], m3[k], c, k);

	for(k = lb+1; k <= lt-1; k++){
		j = k-1;
		zero3(res, u + ir[k] * sizeof(double), m1[k], m2[k], m3[k]);
		interp(res, u + ir[j] * sizeof(double), m1[j], m2[j], m3[j], u + ir[k] * sizeof(double), m1[k], m2[k], m3[k], k);
		// zero3(res, &u[ir[k]], m1[k], m2[k], m3[k]);
		// interp(res, &u[ir[j]], m1[j], m2[j], m3[j], &u[ir[k]], m1[k], m2[k], m3[k], k);
		resid(res, u + ir[k] * sizeof(double), r + ir[k] * sizeof(double), r + ir[k] * sizeof(double), m1[k], m2[k], m3[k], a, k);	
		// resid(res, &u[ir[k]], &r[ir[k]], &r[ir[k]], m1[k], m2[k], m3[k], a, k);	
		psinv(res, r + ir[k] * sizeof(double), u + ir[k] * sizeof(double), m1[k], m2[k], m3[k], c, k);
		//printf("k = %d\n", k);
		// psinv(res, &r[ir[k]], &u[ir[k]], m1[k], m2[k], m3[k], c, k);
	}
	
	j = lt - 1;
	k = lt;
	//interp(res, &u[ir[j]], m1[j], m2[j], m3[j], u, n1, n2, n3, k);	
	interp(res, u + ir[j] * sizeof(double), m1[j], m2[j], m3[j], u, n1, n2, n3, k);	
	resid(res, u, v, r, n1, n2, n3, a, k);	
	psinv(res, r, u, n1, n2, n3, c, k);
	// FILE *file_v, *file_r, *file_u;
    // file_v = fopen("v.txt", "w"); 
	// file_r = fopen("r.txt", "w");
	// file_u = fopen("u.txt", "w");

    // if (file_v == NULL) {
    //     perror("Error opening file");
    // }
	// if (file_u == NULL) {
    //     perror("Error opening file");
    // }
	// if (file_r == NULL) {
    //     perror("Error opening file");
    // }
	
	// int cnt = ceil(NR * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(res, u + i * cpupool_mem_size, res->buf, NR * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(res, u + i * cpupool_mem_size, res->buf, cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NR * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_u, "%f \n", res->buf[j]);
	// 	}
	// }
	
	
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(res, r + i * cpupool_mem_size, res->buf, NR * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(res, r + i * cpupool_mem_size, res->buf, cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NR * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_r, "%f \n", res->buf[j]);
	// 	}
	// }

	// cnt = ceil(NV * sizeof(double) / cpupool_mem_size) + 1;
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(res, v + i * cpupool_mem_size, res->buf, NV * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(res, v + i * cpupool_mem_size, res->buf, cpupool_mem_size);
	// 		if (poll_completion(res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < min(NV * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_v, "%f \n", res->buf[j]);
	// 	}
	// }
	// fclose(file_u);
	// fclose(file_v);
	// fclose(file_r);
	// exit(0);


}

/*
 * ---------------------------------------------------------------------
 * norm2u3 evaluates approximations to the l2 norm and the
 * uniform (or l-infinity or chebyshev) norm, under the
 * assumption that the boundaries are periodic or zero. add the
 * boundaries in with half weight (quarter weight on the edges
 * and eighth weight at the corners) for inhomogeneous boundaries.
 * ---------------------------------------------------------------------
 */
static void norm2u3(struct resources *res, uint64_t r, int n1, int n2, int n3, double* rnm2, double* rnmu, int nx, int ny, int nz){
// #ifdef __clang__
// 		using custom_cast = double (*)[n2][n1];
// 		custom_cast r = reinterpret_cast<custom_cast>(pointer_r);
// #else
// 		double (*r)[n2][n1] = (double (*)[n2][n1])pointer_r;
// #endif	

	static double s, rnmu_local;
	double a;
	int i3, i2, i1;

	double dn;

	if(timeron){
		#pragma omp master
			timer_start(T_NORM2);
	}
	dn = 1.0*nx*ny*nz;

	#pragma omp single
	{
		s = 0.0;
		rnmu_local = 0.0;
	}

	#pragma omp for reduction(+:s,rnmu_local) 
	for(i3 = 1; i3 < n3-1; i3++){
		int tid = omp_get_thread_num(); 
		int num_element = min(n2 * n1, static_cast<int>(block_size / sizeof(double)));
		for (int i = 0; i < ceil(n2 * n1 / num_element); i++) {
			myread(res, r + (i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double),  sizeof(double) * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
			for(i2 = 1; i2 < n2-1; i2++){
				for(i1 = 1; i1 < n1-1; i1++){
					s += res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1] * res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1];
					a = fabs(res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1]);
					// s = s + r[i3][i2][i1] * r[i3][i2][i1];
					// a = fabs(r[i3][i2][i1]);
					if(a > rnmu_local) {
						rnmu_local = a;
					}
				}
			}
		}
		// for(i2 = 1; i2 < n2-1; i2++){
		// 	for(i1 = 1; i1 < n1-1; i1++){
		// 		s = s + r[i3][i2][i1] * r[i3][i2][i1];
		// 		a = fabs(r[i3][i2][i1]);
		// 		if(a > rnmu_local){rnmu_local = a;}
		// 	}
		// }
	}
	*rnmu = rnmu_local;
	*rnm2 = sqrt(s/dn);
	if(timeron){
		#pragma omp master
			timer_stop(T_NORM2);
	}
}

/*
 * ---------------------------------------------------------------------
 * power raises an integer, disguised as a double
 * precision real, to an integer power
 * ---------------------------------------------------------------------
 */
static double power(double a, int n){
	double aj;
	int nj;
	double rdummy;
	double power;

	power = 1.0;
	nj = n;
	aj = a;

	while(nj != 0){
		if((nj%2)==1){rdummy = randlc(&power, aj);}
		rdummy = randlc(&aj, aj);
		nj = nj/2;
	}

	return power;
}

/*
 * --------------------------------------------------------------------
 * psinv applies an approximate inverse as smoother: u = u + Cr
 * 
 * this  implementation costs  15A + 4M per result, where
 * A and M denote the costs of Addition and Multiplication.  
 * presuming coefficient c(3) is zero (the NPB assumes this,
 * but it is thus not a general case), 2A + 1M may be eliminated,
 * resulting in 13A + 3M.
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void psinv(struct resources *res, uint64_t r, uint64_t u, int n1, int n2, int n3, double c[4], int k){
// #ifdef __clang__
// 	using custom_cast = double (*)[n2][n1];
// 	custom_cast r = reinterpret_cast<custom_cast>(pointer_r);	
// 	using custom_cast2 = double (*)[n2][n1];
// 	custom_cast2 u = reinterpret_cast<custom_cast2>(pointer_u);
// #else
// 	double (*r)[n2][n1] = (double (*)[n2][n1])pointer_r;
// 	double (*u)[n2][n1] = (double (*)[n2][n1])pointer_u;	
// #endif		
	int i3, i2, i1;
	double r1[M], r2[M];

	if(timeron){
		#pragma omp master
			timer_start(T_PSINV);
	}
	
	#pragma omp for
	for(i3 = 1; i3 < n3-1; i3++){
		int tid = omp_get_thread_num(); 
		//int num_element = min(4 * n2 * n1, static_cast<int>(block_size / sizeof(double)));
		//for (int i = 0; i < ceil(4 * n2 * n1 / num_element) + 1; i++) {
			myread(res, r + (i3 - 1) * n2 * n1 * sizeof(double), res->buf,  sizeof(double) * 3 * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
			myread(res, u + (i3 * n2 * n1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 3 * n2 * n1,  sizeof(double) * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}

			for(i2 = 1; i2 < n2-1; i2++){
				for(i1 = 0; i1 < n1; i1++){
					r1[i1] = res->buf[n2 * n1 + (i2 - 1) * n1 + i1] + res->buf[n2 * n1 + (i2 + 1) * n1 + i1]
						+ res->buf[i2 * n1 + i1] + res->buf[2 * n2 * n1  + i2 * n1 + i1];
				// r1[i1] = r[i3][i2-1][i1] + r[i3][i2+1][i1]
				// 	+ r[i3-1][i2][i1] + r[i3+1][i2][i1];
				// r2[i1] = r[i3-1][i2-1][i1] + r[i3-1][i2+1][i1]
				// 	+ r[i3+1][i2-1][i1] + r[i3+1][i2+1][i1];
					r2[i1] = res->buf[(i2 - 1) * n1 + i1] + res->buf[(i2 + 1) * n1 + i1]
						+ res->buf[2 * n2 * n1 + (i2 - 1) * n1 + i1] + res->buf[2 * n2 * n1 + (i2 + 1) * n1 + i1];
				}
				for(i1 = 1; i1 < n1-1; i1++){
				
				// u[i3][i2][i1] = u[i3][i2][i1]
				// 	+ c[0] * r[i3][i2][i1]
				// 	+ c[1] * ( r[i3][i2][i1-1] + r[i3][i2][i1+1] + r1[i1] )
				// 	+ c[2] * ( r2[i1] + r1[i1-1] + r1[i1+1] );
					res->buf[3 * n2 * n1 + i2 * n1 + i1] += c[0] * res->buf[n2 * n1 + i2 * n1 + i1]	 	
														  + c[1] * (res->buf[n2 * n1 + i2 * n1 + i1 - 1] + res->buf[n2 * n1 + i2 * n1 + i1 + 1] + r1[i1])
														  + c[2] * (r2[i1] + r1[i1 - 1] + r1[i1 + 1]);
					
					mywrite(res, u + (i3 * n2 * n1 + i2 * n1 + i1) * sizeof(double), res->buf + 3 * n2 * n1 + i2 * n1 + i1,  sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				/*
				 * --------------------------------------------------------------------
				 * assume c(3) = 0    (enable line below if c(3) not= 0)
				 * --------------------------------------------------------------------
				 * > + c(3) * ( r2(i1-1) + r2(i1+1) )
				 * --------------------------------------------------------------------
				 */
				}
			}
			
			// mywrite(res, r + ((i3 - 1) * n2 * n1 + i * n2 * n1) * sizeof(double), res->buf + tid * block_size / sizeof(double),  sizeof(double) * 3 * n2 * n1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
			// mywrite(res, u + (i3 * n2 * n1 + i * n2 * n1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 3 * n2 * n1,  sizeof(double) * n2 * n1);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
	}
	if(timeron){
		#pragma omp master
			timer_stop(T_PSINV);
	}

	/*
	 * --------------------------------------------------------------------
	 * exchange boundary points
	 * --------------------------------------------------------------------
	 */
	//comm3(res, u,n1,n2,n3,k);

	if(debug_vec[0] >= 1){
		#pragma omp single
			rep_nrm(res, u,n1,n2,n3,(char*)"   psinv", k);
	}

	if(debug_vec[3] >= k){
		#pragma omp single
			showall(res, u,n1,n2,n3);
	}

}
/*
 * ---------------------------------------------------------------------
 * report on norm
 * ---------------------------------------------------------------------
 */
static void rep_nrm(struct resources *res, uint64_t u, int n1, int n2, int n3, char* title, int kk){
	double rnm2, rnmu;
	norm2u3(res,u,n1,n2,n3,&rnm2,&rnmu,nx[kk],ny[kk],nz[kk]);
	#pragma omp master
		printf(" Level%2d in %8s: norms =%21.14e%21.14e\n", kk, title, rnm2, rnmu);
}

/*
 * --------------------------------------------------------------------
 * resid computes the residual: r = v - Au
 *
 * this  implementation costs  15A + 4M per result, where
 * A and M denote the costs of addition (or subtraction) and 
 * multiplication, respectively. 
 * presuming coefficient a(1) is zero (the NPB assumes this,
 * but it is thus not a general case), 3A + 1M may be eliminated,
 * resulting in 12A + 3M.
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void resid(struct resources *res, uint64_t u, uint64_t v, uint64_t r, int n1, int n2, int n3, double a[4], int k) {
// #ifdef __clang__
// 	using custom_cast = double (*)[n2][n1];
// 	custom_cast u = reinterpret_cast<custom_cast>(pointer_u);	
// 	using custom_cast2 = double (*)[n2][n1];
// 	custom_cast2 v = reinterpret_cast<custom_cast2>(pointer_v);
// 	using custom_cast3 = double (*)[n2][n1];
// 	custom_cast3 r = reinterpret_cast<custom_cast3>(pointer_r);	
// #else
// 	double (*u)[n2][n1] = (double (*)[n2][n1])pointer_u;
// 	double (*v)[n2][n1] = (double (*)[n2][n1])pointer_v;
// 	double (*r)[n2][n1] = (double (*)[n2][n1])pointer_r;		
// #endif

	int i3, i2, i1;
	double u1[M], u2[M];

	if(timeron){
		#pragma omp master
			timer_start(T_RESID);
	}
	#pragma omp for
	for(i3 = 1; i3 < n3-1; i3++){
		int tid = 0; 
		int num_element = min(4 * n2 * n1, static_cast<int>(block_size / sizeof(double)));
		for (int i = 0; i < 0 + 1; i++) {
			myread(res, u + ((i3 - 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size /sizeof(double),  sizeof(double) * 3 * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
			myread(res, v + (i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 3 * n2 * n1,  sizeof(double) * n2 * n1);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}

			for(i2 = 1; i2 < n2-1; i2++){
				for(i1 = 0; i1 < n1; i1++){
					u1[i1] = res->buf[tid * block_size / sizeof(double) + n2 * n1 + (i2 - 1) * n1 + i1] + res->buf[tid * block_size / sizeof(double) + n2 * n1 + (i2 + 1) * n1 + i1]
							+ res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1] + res->buf[tid * block_size / sizeof(double) + 2 * n2 * n1 + i2 * n1 + i1];
					// u1[i1] = u[i3][i2-1][i1] + u[i3][i2+1][i1]
					// 	+ u[i3-1][i2][i1] + u[i3+1][i2][i1];
					// u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
					// 	+ u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
					u2[i1] = res->buf[tid * block_size / sizeof(double) + (i2 - 1) * n1 + i1] + res->buf[tid * block_size / sizeof(double) + (i2 + 1) * n1 + i1]
							+ res->buf[tid * block_size  / sizeof(double) + 2 * n2 * n1 + (i2 - 1) * n1 + i1] + res->buf[tid * block_size / sizeof(double) + 2 * n2 * n1 + (i2 + 1) * n1 + i1];
				}
			
			for(i1 = 1; i1 < n1-1; i1++){
				/*reuse block*//******************************************************************************************************************************/
				res->buf[tid * block_size / sizeof(double) + 4 * n2 * n1 + i2 * n1 + i1] = res->buf[tid * block_size / sizeof(double) + 3 * n2 * n1 + i2 * n1 + i1] 
																						- a[0] * res->buf[tid * block_size / sizeof(double) + n2 * n1 + i2 * n1 + i1] 
																						- a[2] * (u2[i1] + u1[i1-1] + u1[i1+1]) 
																						- a[3] * (u2[i1-1] + u2[i1 + 1]);
				mywrite(res, r + (i3 * n1 * n2 + i2 * n1 + i1) * sizeof(double), res->buf + tid * block_size / sizeof(double) + 4 * n2 * n1 + i2 * n1 + i1, sizeof(double));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				// r[i3][i2][i1] = v[i3][i2][i1]
				// 	- a[0] * u[i3][i2][i1]
				// 	/*
				// 	 * ---------------------------------------------------------------------
				// 	 * assume a(1) = 0 (enable 2 lines below if a(1) not= 0)
				// 	 * ---------------------------------------------------------------------
				// 	 * > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
				// 	 * > + u1(i1) )
				// 	 * ---------------------------------------------------------------------
				// 	 */
				// 	- a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
				// 	- a[3] * ( u2[i1-1] + u2[i1 + 1]);
				}
			}
		}
		// for(i2 = 1; i2 < n2-1; i2++){
		// 	for(i1 = 0; i1 < n1; i1++){
		// 		u1[i1] = u[i3][i2-1][i1] + u[i3][i2+1][i1]
		// 			+ u[i3-1][i2][i1] + u[i3+1][i2][i1];
		// 		u2[i1] = u[i3-1][i2-1][i1] + u[i3-1][i2+1][i1]
		// 			+ u[i3+1][i2-1][i1] + u[i3+1][i2+1][i1];
		// 	}
		// 	for(i1 = 1; i1 < n1-1; i1++){
		// 		r[i3][i2][i1] = v[i3][i2][i1]
		// 			- a[0] * u[i3][i2][i1]
		// 			/*
		// 			 * ---------------------------------------------------------------------
		// 			 * assume a(1) = 0 (enable 2 lines below if a(1) not= 0)
		// 			 * ---------------------------------------------------------------------
		// 			 * > - a(1) * ( u(i1-1,i2,i3) + u(i1+1,i2,i3)
		// 			 * > + u1(i1) )
		// 			 * ---------------------------------------------------------------------
		// 			 */
		// 			- a[2] * ( u2[i1] + u1[i1-1] + u1[i1+1] )
		// 			- a[3] * ( u2[i1-1] + u2[i1+1] );
		// 	}
		// }
	}
	if(timeron){
		#pragma omp master
			timer_stop(T_RESID);
	}

	/*
	 * --------------------------------------------------------------------
	 * exchange boundary data
	 * --------------------------------------------------------------------
	 */
	//comm3(res, r, n1,n2,n3,k);

	if(debug_vec[0] >= 1){
		#pragma omp single
			rep_nrm(res, r,n1,n2,n3,(char*)"   resid",k);
	}

	if(debug_vec[2] >= k){
		#pragma omp single
			showall(res, r,n1,n2,n3);
	}
}

/*
 * --------------------------------------------------------------------
 * rprj3 projects onto the next coarser grid, 
 * using a trilinear finite element projection: s = r' = P r
 *     
 * this  implementation costs 20A + 4M per result, where
 * A and M denote the costs of addition and multiplication.  
 * note that this vectorizes, and is also fine for cache 
 * based machines.  
 * --------------------------------------------------------------------
 */
static void rprj3(struct resources *res, uint64_t r, int m1k, int m2k, int m3k, uint64_t s, int m1j, int m2j, int m3j, int k){
// #ifdef __clang__
// 	using custom_cast = double (*)[m2k][m1k];
// 	custom_cast r = reinterpret_cast<custom_cast>(pointer_r);
// 	using custom_cast2 = double (*)[m2j][m1j];
// 	custom_cast2 s = reinterpret_cast<custom_cast2>(pointer_s);
// #else
// 	double (*r)[m2k][m1k] = (double (*)[m2k][m1k])pointer_r;
// 	double (*s)[m2j][m1j] = (double (*)[m2j][m1j])pointer_s;		
// #endif	

	int j3, j2, j1, i3, i2, i1, d1, d2, d3, j;

	double x1[M], y1[M], x2, y2;

	if(timeron){
		#pragma omp master
			timer_start(T_RPRJ3);
	}
	if(m1k == 3){
		d1 = 2;
	}else{
		d1 = 1;	
	}
	if(m2k == 3){
		d2 = 2;
	}else{
		d2 = 1;
	}
	if(m3k == 3){
		d3 = 2;
	}else{
		d3 = 1;
	}
	
	#pragma omp for
	for(j3 = 1; j3 < m3j-1; j3++){
		//int tid = omp_get_thread_num(); 
		i3 = 2*j3-d3;	
		//int num_element = min(4 * NM * NM, static_cast<int>(block_size / sizeof(double)));
		
		//for (int i = 0; i < ceil(NM * NM / num_element) + 1; i++) {
			myread(res, r + (i3 * m1k * m2k) * sizeof(double), res->buf,  sizeof(double) * 3 * m1k * m2k);
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
			
			// myread(res, s + (j3 * m1j * m2j ) * sizeof(double), res->buf + 3 * m1k * m2k,  sizeof(double) * m1j * m2j);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }

			for(j2 = 1; j2 < m2j-1; j2++){
				i2 = 2*j2-d2;
				for(j1 = 1; j1 < m1j; j1++){
					i1 = 2 * j1-d1;	

					x1[i1] = res->buf[m1k * m2k + i2 * m1k + i1] + res->buf[m1k * m2k + (i2 + 2) * m1k + i1]	
							+ res->buf[(i2 + 1) * m1k + i1] + res->buf[2 * m1k * m2k + (i2 + 1) * m1k + i1];
					y1[i1] = res->buf[i2 * m1k + i1] + res->buf[2 * m1k * m2k + i2 * m1k + i1]	
							+ res->buf[(i2 + 2) * m1k + i1] + res->buf[2 * m1k * m2k + (i2 + 2) * m1k + i1];
					// x1[i1] = r[i3+1][i2][i1] + r[i3+1][i2+2][i1]
					// 	+ r[i3][i2+1][i1] + r[i3+2][i2+1][i1];

					// y1[i1] = r[i3][i2][i1] + r[i3+2][i2][i1]
					// 	+ r[i3][i2+2][i1] + r[i3+2][i2+2][i1];
				}
				for(j1 = 1; j1 < m1j-1; j1++){
					i1 = 2*j1-d1;
					y2 = res->buf[i2 * m1k + i1 + 1] + res->buf[2 * m1k * m2k + i2 * m1k + i1 + 1]
							+ res->buf[(i2 + 2) * m1k + i1 + 1] + res->buf[2 * m1k * m2k + (i2 + 2) * m1k + i1 + 1];
					x2 = res->buf[m1k * m2k + i2 * m1k + i1 + 1] + res->buf[m1k * m2k  + (i2 + 2) * m1k + i1 + 1]	
							+ res->buf[(i2 + 1) * m1k + i1 + 1] + res->buf[2 * m1k * m2k + (i2 + 1) * m1k + i1 + 1];		
								
					// y2 = r[i3][i2][i1+1] + r[i3+2][i2][i1+1]
					// 	+ r[i3][i2+2][i1+1] + r[i3+2][i2+2][i1+1];

					// x2 = r[i3+1][i2][i1+1] + r[i3+1][i2+2][i1+1]
					// 	+ r[i3][i2+1][i1+1] + r[i3+2][i2+1][i1+1];

					res->buf[3 * m1k * m2k  + j2 * m1j + j1]
					/*s[j3][j2][j1]*/ = 0.5 * res->buf[m1k * m2k + (i2 + 1) * m1k + i1 + 1] 
									+ 0.25 * (res->buf[m1k * m2k  + (i2 + 1) * m1k + i1] + res->buf[m1k * m2k + (i2 + 1) * m1k + i1 + 2] + x2)
									+ 0.125 * (x1[i1] + x1[i1 + 2] + y2)
									+0.0625 * (y1[i1] + y1[i1+2]);
						// 0.5 * r[i3+1][i2+1][i1+1]
						// + 0.25 * ( r[i3+1][i2+1][i1] + r[i3+1][i2+1][i1+2] + x2)
						// + 0.125 * ( x1[i1] + x1[i1+2] + y2)
						// + 0.0625 * ( y1[i1] + y1[i1+2] );
					//printf("%f \n", res->buf[3 * m1k * m2k  + j2 * m1j + j1]);
					mywrite(res, s + (j3 * m1j * m2j + j2 * m1j + j1) * sizeof(double), res->buf + 3 * m1k * m2k + j2 * m1j + j1,  sizeof(double));
					if (poll_completion(res)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}
			}	
			// mywrite(res, r + (i3 * m1k * m2k) * sizeof(double), res->buf,  sizeof(double) * 3 * m1k * m2k);
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
			
	} 
	if(timeron){
		#pragma omp master
			timer_stop(T_RPRJ3);
	}

	j=k-1;



	//comm3(res, s, m1j, m2j, m3j,j);
	// for (int i = 0; i < m3j; i++) {
	// 	myread(res, s + (i * m1j * m2j ) * sizeof(double), res->buf + 3 * m1k * m2k,  sizeof(double) * m1j * m2j);
	// 	if (poll_completion(res)) {
	// 		fprintf(stderr, "poll completion failed\n");
	// 	}
	// 	for (int j = 0; j < m2j;  j++) {
	// 		for (int k = 0; k < m1j; k++) {
	// 			printf("%f \n", res->buf[3 * m1k * m2k + j * m1j + k]);
	// 		}
	// 	}
	// }
	// exit(0);
	
	if(debug_vec[0] >= 1){
		#pragma omp single
			rep_nrm(res,s,m1j,m2j,m3j,(char*)"   rprj3",k-1);	
	}

	if(debug_vec[4] >= k){
		#pragma omp single
			showall(res,s,m1j,m2j,m3j);
	}
}

static void setup(int* n1, int* n2, int* n3, int k){
	int j;

	int ax, mi[MAXLEVEL+1][3];
	int ng[MAXLEVEL+1][3];

	ng[lt][0] = nx[lt];
	ng[lt][1] = ny[lt];
	ng[lt][2] = nz[lt];
	for(ax = 0; ax < 3; ax++){
		for(k = lt-1; k >= 1; k--){
			ng[k][ax] = ng[k+1][ax]/2;
		}
	}
	for(k = lt; k >= 1; k--){
		nx[k] = ng[k][0];
		ny[k] = ng[k][1];
		nz[k] = ng[k][2];
	}

	for(k = lt; k >= 1; k--){
		for (ax = 0; ax < 3; ax++){
			mi[k][ax] = 2 + ng[k][ax];
		}

		m1[k] = mi[k][0];
		m2[k] = mi[k][1];
		m3[k] = mi[k][2];
	}

	k = lt;
	is1 = 2 + ng[k][0] - ng[lt][0];
	ie1 = 1 + ng[k][0];
	*n1 = 3 + ie1 - is1;
	is2 = 2 + ng[k][1] - ng[lt][1];
	ie2 = 1 + ng[k][1];
	*n2 = 3 + ie2 - is2;
	is3 = 2 + ng[k][2] - ng[lt][2];
	ie3 = 1 + ng[k][2];
	*n3 = 3 + ie3 - is3;

	ir[lt] = 0;
	for(j = lt-1; j >= 1; j--){
		ir[j] = ir[j+1]+ONE*m1[j+1]*m2[j+1]*m3[j+1];
	}

	if(debug_vec[1] >= 1){
		printf(" in setup, \n");
		printf("   k  lt  nx  ny  nz  n1  n2  n3 is1 is2 is3 ie1 ie2 ie3\n");
		printf("%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d%4d\n", 
				k,lt,ng[k][0],ng[k][1],ng[k][2],*n1,*n2,*n3,is1,is2,is3,ie1,ie2,ie3);
	}
}

static void showall(struct resources *res, uint64_t z, int n1, int n2, int n3){
// #ifdef __clang__
// 	using custom_cast = double (*)[n2][n1];
// 	custom_cast z = reinterpret_cast<custom_cast>(pointer_z);	
// #else
// 	double (*z)[n2][n1] = (double (*)[n2][n1])pointer_z;	
// #endif

	int i1,i2,i3;
	int m1, m2, m3;

	m1 = min(n1,18);
	m2 = min(n2,14);
	m3 = min(n3,18);

	printf("\n");

	for(i3 = 0; i3 < m3; i3++){
		mywrite(res, z + i3 * sizeof(m2 * m1) * sizeof(double), res->buf, sizeof(m2 * m1) * sizeof(double));
		if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
		}
		for(i2 = 0; i2 < m2; i2++){
			for(i1 = 0; i1 < m1; i1++){			
				printf("%6.3f", res->buf[i2 * m1 + i1]);
			}
			printf("\n");
		}
		printf(" - - - - - - - \n");
	}
	printf("\n");
}

static void zero3(struct resources *res, uint64_t z, int n1, int n2, int n3){
// #ifdef __clang__
// 		using custom_cast = double (*)[n2][n1];
// 		custom_cast z = reinterpret_cast<custom_cast>(pointer_z);
// #else
// 		double (*z)[n2][n1] = (double (*)[n2][n1])pointer_z;
// #endif

	int i1, i2, i3;
	memset(res->buf, 0.0, n2 * n1 * sizeof(double));

    #pragma omp for  
	for(i3 = 0; i3 < n3; i3++){
		int tid = omp_get_thread_num();
		mywrite(res, z + i3 * n2 * n1 * sizeof(double), res->buf + tid * block_size / sizeof(double), n2 * n1 * sizeof(double));
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
		// for(i2 = 0; i2 < n2; i2++){
		// 	for(i1 = 0; i1 < n1; i1++){
		// 		z[i3][i2][i1] = 0.0;
		// 	}
		// }
	}
}

/*
 * ---------------------------------------------------------------------
 * zran3 loads +1 at ten randomly chosen points,
 * loads -1 at a different ten random points,
 * and zero elsewhere.
 * ---------------------------------------------------------------------
 */
static void zran3(struct resources *res, uint64_t z, int n1, int n2, int n3, int nx, int ny, int k){
// #ifdef __clang__
// 	using custom_cast = double (*)[n2][n1];
// 	custom_cast z = reinterpret_cast<custom_cast>(pointer_z);
// #else
// 	double (*z)[n2][n1] = (double (*)[n2][n1])pointer_z;
// #endif

	int i0, m0, m1;

	int i1, i2, i3, d1, e1, e2, e3;
	double xx, x0, x1, a1, a2, ai;

	double ten[2][MM], best;
	int i, j1[2][MM], j2[2][MM], j3[2][MM];
	int jg[2][MM][4];

	a1 = power(A, nx);
	a2 = power(A, nx*ny);

	#pragma omp parallel
		zero3(res, z, n1, n2, n3);

	

	i = is1-2+nx*(is2-2+ny*(is3-2));

	ai = power(A, i); 
	d1 = ie1 - is1 + 1;
	e1 = ie1 - is1 + 2;
	e2 = ie2 - is2 + 2;
	e3 = ie3 - is3 + 2;
	x0 = X;
	randlc(&x0, ai);
	//double* ptr = (double*)(uintptr_t)z;
	for(i3 = 1; i3 < e3; i3++){
		x1 = x0;
		for(i2 = 1; i2 < e2; i2++){
			xx = x1;
			vranlc(res, d1, &xx, A, z + (i3 * n2 * n1 + i2 * n1 + 1) * sizeof(double)/*&(z[i3][i2][1])*/);
			randlc(&x1,a1);
		}
		randlc(&x0, a2);
	}


	
	
	/*
	 * ---------------------------------------------------------------------
	 * each processor looks for twenty candidates
	 * ---------------------------------------------------------------------
	 */	
	for(i = 0; i < MM; i++){
		ten[1][i] = 0.0;
		j1[1][i] = 0;
		j2[1][i] = 0;
		j3[1][i] = 0;
		ten[0][i] = 1.0;
		j1[0][i] = 0;
		j2[0][i] = 0;
		j3[0][i] = 0;
	}

	int tid = omp_get_thread_num(); 
	for(i3 = 1; i3 < n3-1; i3++){
		myread(res, z + i3 * n2 * n1 * sizeof(double), res->buf + tid * block_size / sizeof(double), sizeof(double) * n2 * n1);
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
		for(i2 = 1; i2 < n2-1; i2++){
			for(i1 = 1; i1 < n1-1; i1++){
				//printf("z: %f\n", res->buf[tid * block_size + i2 * n1 + i1]);
				if(/*z[i3][i2][i1]*/ res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1] > ten[1][0]){
					ten[1][0] = res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1]/*z[i3][i2][i1]*/;
					j1[1][0] = i1;
					j2[1][0] = i2;
					j3[1][0] = i3;
					bubble(ten, j1, j2, j3, MM, 1);
				}
				if(/*z[i3][i2][i1]*/ res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1] < ten[0][0]){
					ten[0][0] = res->buf[tid * block_size / sizeof(double) + i2 * n1 + i1]/*z[i3][i2][i1]*/;
					j1[0][0] = i1;
					j2[0][0] = i2;
					j3[0][0] = i3;
					bubble(ten, j1, j2, j3, MM, 0);
				}
			}
		}
	}

	/*
	 * ---------------------------------------------------------------------
	 * now which of these are globally best?
	 * ---------------------------------------------------------------------
	 */	
	i1 = MM - 1;
	i0 = MM - 1; 
	for(i = MM - 1; i >= 0; i--){
		best = 0.0;
		if(best < ten[1][i1]){
			jg[1][i][0] = 0;
			jg[1][i][1] = is1 - 2 + j1[1][i1];
			jg[1][i][2] = is2 - 2 + j2[1][i1];
			jg[1][i][3] = is3 - 2 + j3[1][i1];
			i1 = i1-1;
		}else{
			jg[1][i][0] = 0;
			jg[1][i][1] = 0;
			jg[1][i][2] = 0;
			jg[1][i][3] = 0;
		}
		best = 1.0;
		if(best > ten[0][i0]){
			jg[0][i][0] = 0;
			jg[0][i][1] = is1 - 2 + j1[0][i0];
			jg[0][i][2] = is2 - 2 + j2[0][i0];
			jg[0][i][3] = is3 - 2 + j3[0][i0];
			i0 = i0-1;
		}else{
			jg[0][i][0] = 0;
			jg[0][i][1] = 0;
			jg[0][i][2] = 0;
			jg[0][i][3] = 0;
		}
	}
	m1 = 0;
	m0 = 0;
	memset(res->buf, 0.0, n2 * n1 * sizeof(double));

	for(i3 = 0; i3 < n3; i3++){
		mywrite(res, z + i3 * n2 * n1 * sizeof(double), res->buf + tid * block_size / sizeof(double), n2 * n1 * sizeof(double));
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
		// for(i2 = 0; i2 < n2; i2++){
		// 	for(i1 = 0; i1 < n1; i1++){
		// 		z[i3][i2][i1] = 0.0;
		// 	}
		// }
	}

	//printf("z = %ld\n", z);
	for (i = MM-1; i >= m0; i--){
		//z[jg[0][i][3]][jg[0][i][2]][jg[0][i][1]] = -1.0;
		res->buf[0] = -1.0;
		//printf("jg -1 = %d, %d, %d\n", jg[0][i][3], jg[0][i][2], jg[0][i][1]);
		//printf("z addr is = %ld\n", z + (jg[0][i][3] * n2 * n1 + jg[0][i][2] * n1 + jg[0][i][1]) * sizeof(double));
		//printf("jg = %d\n", jg[0][i][3] * n2 * n1 + jg[0][i][2] * n1 + jg[0][i][1]);
		mywrite(res, z + (jg[0][i][3] * n2 * n1 + jg[0][i][2] * n1 + jg[0][i][1]) * sizeof(double), &res->buf[0], sizeof(double));
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n\n");
		}
	} 
	//printf("z1 = %ld\n", z);

	printf("\n\n");

	for(i = MM-1; i >= m1; i--){
		//printf("jg +1 = %d, %d, %d\n", jg[1][i][3], jg[1][i][2], jg[1][i][1]);
		//z[jg[1][i][3]][jg[1][i][2]][jg[1][i][1]] = +1.0;
		res->buf[tid * block_size / sizeof(double)] = +1.0;
		mywrite(res, z + (jg[1][i][3] * n2 * n1 + jg[1][i][2] * n1 + jg[1][i][1]) * sizeof(double), res->buf + tid * block_size / sizeof(double) , sizeof(double));
		if (poll_completion(res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	}

	
	

	//#pragma omp parallel 
	//comm3(res, z, n1, n2, n3, k);
	
	
}

static int sock_connect(const char *servername, int port)
{
	struct addrinfo *resolved_addr = NULL;
	struct addrinfo *iterator;
	char service[6];
	int sockfd = -1;
	int listenfd = 0;
	int tmp;
	struct addrinfo hints =
		{
			.ai_flags = AI_PASSIVE,
			.ai_family = AF_INET,
			.ai_socktype = SOCK_STREAM};
	if (sprintf(service, "%d", port) < 0)
		goto sock_connect_exit;
	/* Resolve DNS address, use sockfd as temp storage */
	sockfd = getaddrinfo(servername, service, &hints, &resolved_addr);
	if (sockfd < 0)
	{
		fprintf(stderr, "%s for %s:%d\n", gai_strerror(sockfd), servername, port);
		goto sock_connect_exit;
	}
	/* Search through results and find the one we want */
	for (iterator = resolved_addr; iterator; iterator = iterator->ai_next)
	{
		sockfd = socket(iterator->ai_family, iterator->ai_socktype, iterator->ai_protocol);
		if (sockfd >= 0)
		{
			if (servername){
				/* Client mode. Initiate connection to remote */
				if ((tmp = connect(sockfd, iterator->ai_addr, iterator->ai_addrlen)))
				{
					fprintf(stdout, "failed connect \n");
					perror("connect failed");
					//exit(0);
					close(sockfd);
					sockfd = -1;
				}
            }
			else
			{
					/* Server mode. Set up listening socket an accept a connection */
					listenfd = sockfd;
					sockfd = -1;
					if (bind(listenfd, iterator->ai_addr, iterator->ai_addrlen))
						goto sock_connect_exit;
					listen(listenfd, 1);
					sockfd = accept(listenfd, NULL, 0);
			}
		}
	}
sock_connect_exit:
	if (listenfd)
		close(listenfd);
	if (resolved_addr)
		freeaddrinfo(resolved_addr);
	if (sockfd < 0)
	{
		if (servername)
			fprintf(stderr, "Couldn't connect to %s:%d\n", servername, port);
		else
		{
			perror("server accept");
			fprintf(stderr, "accept() failed\n");
		}
	}
	return sockfd;
}
/******************************************************************************
* Function: sock_sync_data
*
* Input
* sock socket to transfer data on
* xfer_size size of data to transfer
* local_data pointer to data to be sent to remote
*
* Output
* remote_data pointer to buffer to receive remote data
*
* Returns
* 0 on success, negative error code on failure
*
* Description
* Sync data across a socket. The indicated local data will be sent to the
* remote. It will then wait for the remote to send its data back. It is
* assumed that the two sides are in sync and call this function in the proper
* order. Chaos will ensue if they are not. :)
*
* Also note this is a blocking function and will wait for the full data to be
* received from the remote.
*
******************************************************************************/
int sock_sync_data(int sock, int xfer_size, char *local_data, char *remote_data)
{
	int rc;
	int read_bytes = 0;
	int total_read_bytes = 0;
	rc = write(sock, local_data, xfer_size);
	if (rc < xfer_size)
		fprintf(stderr, "Failed writing data during sock_sync_data\n");
	else
		rc = 0;
	while (!rc && total_read_bytes < xfer_size)
	{
		read_bytes = read(sock, remote_data, xfer_size);
		if (read_bytes > 0)
			total_read_bytes += read_bytes;
		else
			rc = read_bytes;
	}
	return rc;
}
/******************************************************************************
End of socket operations
******************************************************************************/
/* poll_completion */
/******************************************************************************
* Function: poll_completion
*
* Input
* res pointer to resources structure
*
* Output
* none
*
* Returns
* 0 on success, 1 on failure
*
* Description
* Poll the completion queue for a single event. This function will continue to
* poll the queue until MAX_POLL_CQ_TIMEOUT milliseconds have passed.
*
******************************************************************************/
static int poll_completion(struct resources *res)
{
	struct ibv_wc wc[1];
	unsigned long start_time_msec;
	unsigned long cur_time_msec;
	struct timeval cur_time;
	int poll_result = 0;
	int rc = 0;
	/* poll the completion for a while before giving up of doing it .. */
	gettimeofday(&cur_time, NULL);
	start_time_msec = (cur_time.tv_sec * 1000) + (cur_time.tv_usec / 1000);
	do
	{
		poll_result = ibv_poll_cq(res->cq, 1, wc);
		gettimeofday(&cur_time, NULL);
		cur_time_msec = (cur_time.tv_sec * 1000) + (cur_time.tv_usec / 1000);
	} while ((poll_result == 0) && ((cur_time_msec - start_time_msec) < MAX_POLL_CQ_TIMEOUT));
	if (poll_result < 0)
	{
		/* poll CQ failed */
		fprintf(stderr, "poll CQ failed\n");
		rc = 1;
	}
	else if (poll_result == 0)
	{ /* the CQ is empty */
		fprintf(stderr, "completion wasn't found in the CQ after timeout\n");
		rc = 1;
	}
	else
	{
		/* CQE found */
		//fprintf(stdout, "completion was found in CQ with status 0x%x\n", wc.status);
		for (int i = 0; i < 1; i++) {
			if (wc[i].status != IBV_WC_SUCCESS)
			{
				printf("got bad completion with status: 0x%x, vendor syndrome: 0x%x\n", wc[i].status, wc[i].vendor_err);
				rc = 1;
			}
		}
		
	}
	return rc;
}

static int poll_completion_2(struct resources *res)
{
	struct ibv_wc wc;
	unsigned long start_time_msec;
	unsigned long cur_time_msec;
	struct timeval cur_time;
	int poll_result;
	int rc = 0;
	/* poll the completion for a while before giving up of doing it .. */
	gettimeofday(&cur_time, NULL);
	start_time_msec = (cur_time.tv_sec * 1000) + (cur_time.tv_usec / 1000);
	do
	{
		poll_result = ibv_poll_cq(res->cq, 2, &wc);
		gettimeofday(&cur_time, NULL);
		cur_time_msec = (cur_time.tv_sec * 1000) + (cur_time.tv_usec / 1000);
	} while ((poll_result == 0) && ((cur_time_msec - start_time_msec) < MAX_POLL_CQ_TIMEOUT));
	if (poll_result < 0)
	{
		/* poll CQ failed */
		fprintf(stderr, "poll CQ failed\n");
		rc = 1;
	}
	else if (poll_result == 0)
	{ /* the CQ is empty */
		fprintf(stderr, "completion wasn't found in the CQ after timeout\n");
		rc = 1;
	}
	else
	{
		/* CQE found */
		//fprintf(stdout, "completion was found in CQ with status 0x%x\n", wc.status);
		/* check the completion status (here we don't care about the completion opcode */
		if (wc.status != IBV_WC_SUCCESS)
		{
			fprintf(stderr, "got bad completion with status: 0x%x, vendor syndrome: 0x%x\n", wc.status, wc.vendor_err);
			exit(0);
			rc = 1;
		}
	}
	return rc;
}


// static int post_malloc(struct resources *res, size_t size)
// {
// 	struct ibv_send_wr sr;
// 	struct ibv_sge sge;
// 	struct ibv_send_wr *bad_wr = NULL;
// 	int rc;
// 	/* prepare the scatter/gather entry */
// 	memset(&sge, 0, sizeof(sge));
// 	sge.addr = (uintptr_t)res->buf;
// 	sge.length = block_size;
// 	sge.lkey = res->mr->lkey;
// 	/* prepare the send work request */
// 	memset(&sr, 0, sizeof(sr));
// 	sr.next = NULL;
// 	sr.wr_id = 0;
// 	sr.sg_list = &sge;
// 	sr.num_sge = 1;
// 	sr.opcode = IBV_WR_SEND;
// 	sr.send_flags = IBV_SEND_SIGNALED;
// 	if (opcode == IBV_WR_ATOMIC_CMP_AND_SWP) {
//         // 原子比较并交换操作
//         sr.wr.atomic.remote_addr = res->remote_props.addr;
//         sr.wr.atomic.rkey = res->remote_props.rkey;
//         sr.wr.atomic.compare_add = res->atomic_compare;
//         sr.wr.atomic.swap = res->atomic_swap;
//     } else if (opcode == IBV_WR_ATOMIC_FETCH_AND_ADD) {
//         // 原子取并加操作
//         sr.wr.atomic.remote_addr = res->remote_props.addr;
//         sr.wr.atomic.rkey = res->remote_props.rkey;
//         sr.wr.atomic.compare_add = res->atomic_add;
// 	} else if (opcode != IBV_WR_SEND) {
// 		sr.wr.rdma.remote_addr = res->remote_props.addr;
//         sr.wr.rdma.rkey = res->remote_props.rkey;
// 	}

// 	/* there is a Receive Request in the responder side, so we won't get any into RNR flow */
// 	rc = ibv_post_send(res->qp, &sr, &bad_wr);
// 	if (rc) {
// 		fprintf(stderr, "failed1 to post SR\n");
// 	}
//     else {
// 		printf(stdout, "Send Request was posted\n");
// 	}
// 		return rc;
// 	}
	



/******************************************************************************
* Function: mymalloc
******************************************************************************/

uint64_t mymalloc(struct resources *res, size_t mem_size, size_t offset) {
	struct ibv_send_wr sr;
	struct ibv_sge sge;
	struct ibv_send_wr *bad_wr = NULL;
	int rc;

	/* prepare the scatter/gather entry */
	memset(&sge, 0, sizeof(sge));
	sge.addr = (uintptr_t)res->buf;
	sge.length = mem_size;
	sge.lkey = res->mr->lkey;
	/* prepare the send work request */
	memset(&sr, 0, sizeof(sr));
	sr.next = NULL;
	sr.wr_id = 0;
	sr.sg_list = &sge;
	sr.num_sge = 1;
	sr.opcode = IBV_WR_RDMA_WRITE;
	sr.send_flags = IBV_SEND_SIGNALED; 
	memset(res->buf, 0.0, mem_size);
	
	//set the remote addresss
	sr.wr.rdma.remote_addr = res->remote_props.addr + offset;
    sr.wr.rdma.rkey = res->remote_props.rkey;
	//fprintf(stdout, "Remote Address = %ld" "\n", sr.wr.rdma.remote_addr);
	//fprintf(stdout, "Remote Address = 0x%" PRIx64 "\n", sr.wr.rdma.remote_addr);
	/* there is a Receive Request in the responder side, so we won't get any into RNR flow */
	rc = ibv_post_send(res->qp, &sr, &bad_wr);
	if (rc) {
		fprintf(stderr, "failed1 to post Allocation request\n");
	} else {
		//fprintf(stdout, "Allocation Request was posted\n");
	}
	return sr.wr.rdma.remote_addr;
}

/******************************************************************************
* Function: myread
******************************************************************************/
int myread(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size) {
	struct ibv_send_wr sr;
	struct ibv_sge sge;
	struct ibv_send_wr *bad_wr = NULL;
	int rc;

	/* prepare the scatter/gather entry */
	memset(&sge, 0, sizeof(sge));
	//sge.addr = (uintptr_t)res->buf + block_size * index;
	sge.addr = (uintptr_t)local_buf;
	sge.length = size;
	sge.lkey = res->mr->lkey;
	/* prepare the send work request */
	memset(&sr, 0, sizeof(sr));
	sr.next = NULL;
	sr.wr_id = 0;
	sr.sg_list = &sge;
	sr.num_sge = 1;
	sr.opcode = IBV_WR_RDMA_READ;
	sr.send_flags = IBV_SEND_SIGNALED;

	sr.wr.rdma.remote_addr = remote_mem_addr;
    sr.wr.rdma.rkey = res->remote_props.rkey;
	//fprintf(stdout, "Remote Address that will be read = %ld" "\n", sr.wr.rdma.remote_addr);
	//fprintf(stdout, "Remote Address that will be read = 0x%" PRIx64 "\n", sr.wr.rdma.remote_addr);

	rc = ibv_post_send(res->qp, &sr, &bad_wr);
	if (rc) {
		fprintf(stderr, "failed1 to post SR for read and rc = %d\n", rc);
	}
   // else {
	//	fprintf(stderr, "myread success\n");
	//}
		
	return rc;
}

/******************************************************************************
* Function: mywrite
******************************************************************************/
int mywrite(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size) {
	struct ibv_send_wr sr;
	struct ibv_sge sge;
	struct ibv_send_wr *bad_wr = NULL;
	int rc;

	/* prepare the scatter/gather entry */
	memset(&sge, 0, sizeof(sge));
	//sge.addr = (uintptr_t)res->buf + block_size * index;
	sge.addr = (uintptr_t)local_buf;
	sge.length = size;
	sge.lkey = res->mr->lkey;
	/* prepare the send work request */
	memset(&sr, 0, sizeof(sr));
	sr.next = NULL;
	sr.wr_id = 0;
	sr.sg_list = &sge;
	sr.num_sge = 1;
	sr.opcode = IBV_WR_RDMA_WRITE;
	sr.send_flags = IBV_SEND_SIGNALED;

	sr.wr.rdma.remote_addr = remote_mem_addr;
    sr.wr.rdma.rkey = res->remote_props.rkey;

	rc = ibv_post_send(res->qp, &sr, &bad_wr);
	//printf("rc = %d\n", rc);
	if (rc) {
		fprintf(stderr, "failed1 to post SR for write and rc = %d\n",rc);
	}
	// } else {
	// 	fprintf(stderr, "mywrite success\n");
	// }
	return rc;
}





/******************************************************************************
* Function: post_send
*
* Input
* res pointer to resources structure
* opcode IBV_WR_SEND, IBV_WR_RDMA_READ or IBV_WR_RDMA_WRITE
*
* Output
* none
*
* Returns
* 0 on success, error code on failure
*
* Description
* This function will create and post a send work request
******************************************************************************/
// static int post_send(struct resources *res, int opcode)//, int //index)
// {
// 	struct ibv_send_wr sr;
// 	struct ibv_sge sge;
// 	struct ibv_send_wr *bad_wr = NULL;
// 	int rc;

// 	// size_t* indexes = (size_t*)malloc(block_num * sizeof(size_t));
//     // for (size_t i = 0; i < block_num; i++) {
//     //     indexes[i] = i;
//     // }
// 	// shuffle(indexes, block_num);

// 	/* prepare the scatter/gather entry */
// 	memset(&sge, 0, sizeof(sge));
// 	//sge.addr = (uintptr_t)res->buf + block_size * index;
// 	sge.addr = (uintptr_t)res->buf;
// 	sge.length = block_size;
// 	sge.lkey = res->mr->lkey;
// 	/* prepare the send work request */
// 	memset(&sr, 0, sizeof(sr));
// 	sr.next = NULL;
// 	sr.wr_id = 0;
// 	sr.sg_list = &sge;
// 	sr.num_sge = 1;
// 	sr.opcode = opcode;
// 	sr.send_flags = IBV_SEND_SIGNALED;
// 	if (opcode == IBV_WR_ATOMIC_CMP_AND_SWP) {
//         // 原子比较并交换操作
//         sr.wr.atomic.remote_addr = res->remote_props.addr;
//         sr.wr.atomic.rkey = res->remote_props.rkey;
//         sr.wr.atomic.compare_add = res->atomic_compare;
//         sr.wr.atomic.swap = res->atomic_swap;
//     } else if (opcode == IBV_WR_ATOMIC_FETCH_AND_ADD) {
//         // 原子取并加操作
//         sr.wr.atomic.remote_addr = res->remote_props.addr;
//         sr.wr.atomic.rkey = res->remote_props.rkey;
//         sr.wr.atomic.compare_add = res->atomic_add;
// 	} else if (opcode != IBV_WR_SEND) {
//     //     //sr.wr.rdma.remote_addr = res->remote_props.addr + block_size * index;
// 		sr.wr.rdma.remote_addr = res->remote_props.addr;
//         sr.wr.rdma.rkey = res->remote_props.rkey;
// 	}

//         //sr.wr.rdma.remote_addr = res->remote_props.addr + block_size * index;
// 	// sr.wr.rdma.remote_addr = res->remote_props.addr;
//     // sr.wr.rdma.rkey = res->remote_props.rkey;
// 	// fprintf(stdout, "Remote Address = %ld" "\n", res->remote_props.addr);
// 	// fprintf(stdout, "Remote Address = 0x%" PRIx64 "\n", res->remote_props.addr);


// 	//memset(sr.wr.rdma.remote_addr, 23, block_size);
// 	// if (opcode != IBV_WR_SEND)
// 	// {
// 	// 	sr.wr.rdma.remote_addr = res->remote_props.addr;
// 	// 	sr.wr.rdma.rkey = res->remote_props.rkey;
// 	// }
// 	/* there is a Receive Request in the responder side, so we won't get any into RNR flow */
// 	rc = ibv_post_send(res->qp, &sr, &bad_wr);
// 	if (rc) {
// 		fprintf(stderr, "failed1 to post SR\n");
// 	}
//     else {
// 		switch (opcode)
// 		{
// 		case IBV_WR_SEND:
// 			fprintf(stdout, "Send Request was posted\n");
// 			break;
// 		case IBV_WR_RDMA_READ:
// 			//fprintf(stdout, "RDMA Read Request was posted\n");
// 			break;
// 		case IBV_WR_RDMA_WRITE:
// 			//fprintf(stdout, "RDMA Write Request was posted\n");
// 			break;
// 		default:
// 			//fprintf(stdout, "Unknown Request was posted\n");
// 			break;
// 		}
// 	}
// 	return rc;
// }
/******************************************************************************
* Function: post_receive
*
* Input
* res pointer to resources structure
*
* Output
* none
*
* Returns
* 0 on success, error code on failure
*
* Description
*
******************************************************************************/
static int post_receive(struct resources *res)
{
	struct ibv_recv_wr rr;
	struct ibv_sge sge;
	struct ibv_recv_wr *bad_wr;
	int rc;
	/* prepare the scatter/gather entry */
	memset(&sge, 0, sizeof(sge));
	sge.addr = (uintptr_t)res->buf;
	sge.length = 1024;
	sge.lkey = res->mr->lkey;
	/* prepare the receive work request */
	memset(&rr, 0, sizeof(rr));
	rr.next = NULL;
	rr.wr_id = 0;
	rr.sg_list = &sge;
	rr.num_sge = 1;
	/* post the Receive Request to the RQ */
	rc = ibv_post_recv(res->qp, &rr, &bad_wr);
	if (rc)
		fprintf(stderr, "failed to post RR\n");
	else
		fprintf(stdout, "Receive Request was posted\n");
	return rc;
}
/******************************************************************************
* Function: resources_init
*
* Input
* res pointer to resources structure
*
* Output
* res is initialized
*
* Returns
* none
*
* Description
* res is initialized to default values
******************************************************************************/
static void resources_init(struct resources *res)
{
	memset(res, 0, sizeof *res);
	res->sock = -1;
}
/******************************************************************************
* Function: resources_create
*
* Input
* res pointer to resources structure to be filled in
*
* Output
* res filled in with resources
*
* Returns
* 0 on success, 1 on failure
*
* Description
*
* This function creates and allocates all necessary system resources. These
* are stored in res.
*****************************************************************************/
static int resources_create(struct resources *res)
{
	struct ibv_device **dev_list = NULL;
	struct ibv_qp_init_attr qp_init_attr;
	struct ibv_device *ib_dev = NULL;
	//struct ibv_device_attr device_attr;
	size_t size;
	int i;
	int mr_flags = 0;
	int cq_size = 1;
	int num_devices;
	int rc = 0;
	/* if client side */
	if (config.server_name)
	{
		res->sock = sock_connect(config.server_name, config.tcp_port);
		if (res->sock < 0)
		{
			fprintf(stderr, "failed to establish TCP connection to server %s, port %d\n",
					config.server_name, config.tcp_port);
			rc = -1;
			goto resources_create_exit;
		}
	}
	else
	{
		fprintf(stdout, "waiting on port %d for TCP connection\n", config.tcp_port);
		res->sock = sock_connect(NULL, config.tcp_port);
		if (res->sock < 0)
		{
			fprintf(stderr, "failed to establish TCP connection with client on port %d\n",
					config.tcp_port);
			rc = -1;
			goto resources_create_exit;
		}
	}
	fprintf(stdout, "TCP connection was established\n");
	fprintf(stdout, "searching for IB devices in host\n");
	/* get device names in the system */
	dev_list = ibv_get_device_list(&num_devices);
	if (!dev_list)
	{
		fprintf(stderr, "failed to get IB devices list\n");
		rc = 1;
		goto resources_create_exit;
	}
	/* if there isn't any IB device in host */
	if (!num_devices)
	{
		fprintf(stderr, "found %d device(s)\n", num_devices);
		rc = 1;
		goto resources_create_exit;
	}
	fprintf(stdout, "found %d device(s)\n", num_devices);
	/* search for the specific device we want to work with */
	for (i = 0; i < num_devices; i++)
	{
		if (!config.dev_name)
		{
			config.dev_name = strdup(ibv_get_device_name(dev_list[i]));
			fprintf(stdout, "device not specified, using first one found: %s\n", config.dev_name);
		}
		if (!strcmp(ibv_get_device_name(dev_list[i]), config.dev_name))
		{
			ib_dev = dev_list[i];
			break;
		}
	}
	/* if the device wasn't found in host */
	if (!ib_dev)
	{
		fprintf(stderr, "IB device %s wasn't found\n", config.dev_name);
		rc = 1;
		goto resources_create_exit;
	}
	/* get device handle */
	res->ib_ctx = ibv_open_device(ib_dev);
	if (!res->ib_ctx)
	{
		fprintf(stderr, "failed to open device %s\n", config.dev_name);
		rc = 1;
		goto resources_create_exit;
	}
	/* We are now done with device list, free it */
	ibv_free_device_list(dev_list);
	dev_list = NULL;
	ib_dev = NULL;
	/* query port properties */
	if (ibv_query_port(res->ib_ctx, config.ib_port, &res->port_attr))
	{
		fprintf(stderr, "ibv_query_port on port %u failed\n", config.ib_port);
		rc = 1;
		goto resources_create_exit;
	}

	//ibv_query_device(res->ib_ctx, &device_attr);
	//printf("max message size is: %d",res->port_attr.max_msg_sz);
	
	/* allocate Protection Domain */
	res->pd = ibv_alloc_pd(res->ib_ctx);
	if (!res->pd)
	{
		fprintf(stderr, "ibv_alloc_pd failed\n");
		rc = 1;
		goto resources_create_exit;
	}

	/* each side will send only one WR, so Completion Queue with 1 entry is enough */
	cq_size = 500;
	res->cq = ibv_create_cq(res->ib_ctx, cq_size, NULL, NULL, 0);

	if (!res->cq)
	{
		fprintf(stderr, "failed to create CQ with %u entries\n", cq_size);
		rc = 1;
		goto resources_create_exit;
	}

	/* allocate the memory buffer that will hold the data */
	size = cpupool_mem_size;
	res->buf = (double *)malloc(size);
	if (!res->buf)
	{
		fprintf(stderr, "failed to malloc %Zu bytes to memory buffer\n", size);
		rc = 1;
		goto resources_create_exit;
	}
	memset(res->buf, 0.0, size);
	/* only in the client side put the message in the memory buffer */
	// if (config.server_name)
	// {
	// 	strcpy(res->buf, MSG);
	// 	fprintf(stdout, "going to send the message: '%s'\n", res->buf);
	// }
	// else
	// 	memset(res->buf, 0, size);
	// /* register the memory buffer */
	mr_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_ATOMIC;
	//mr_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE;
	res->mr = ibv_reg_mr(res->pd, res->buf, size, mr_flags);
	if (!res->mr)
	{
		fprintf(stderr, "ibv_reg_mr failed with mr_flags=0x%x\n", mr_flags);
		rc = 1;
		goto resources_create_exit;
	}
	fprintf(stdout, "MR was registered with addr=%p, lkey=0x%x, rkey=0x%x, flags=0x%x\n",
			res->buf, res->mr->lkey, res->mr->rkey, mr_flags);
	/* create the Queue Pair */
	memset(&qp_init_attr, 0, sizeof(qp_init_attr));
	qp_init_attr.qp_type = IBV_QPT_RC;
	qp_init_attr.sq_sig_all = 1;
	qp_init_attr.send_cq = res->cq;
	qp_init_attr.recv_cq = res->cq;
	qp_init_attr.cap.max_send_wr = 100;
	qp_init_attr.cap.max_recv_wr = 100;
	qp_init_attr.cap.max_send_sge = 15;
	qp_init_attr.cap.max_recv_sge = 15;
	res->qp = ibv_create_qp(res->pd, &qp_init_attr);
	if (!res->qp)
	{
		fprintf(stderr, "failed to create QP\n");
		rc = 1;
		goto resources_create_exit;
	}
	fprintf(stdout, "QP was created, QP number=0x%x\n", res->qp->qp_num);
resources_create_exit:
	if (rc)
	{
		/* Error encountered, cleanup */
		if (res->qp)
		{
			ibv_destroy_qp(res->qp);
			res->qp = NULL;
		}
		if (res->mr)
		{
			ibv_dereg_mr(res->mr);
			res->mr = NULL;
		}
		if (res->buf)
		{
			free(res->buf);
			res->buf = NULL;
		}
		if (res->cq)
		{
			ibv_destroy_cq(res->cq);
			res->cq = NULL;
		}
		if (res->pd)
		{
			ibv_dealloc_pd(res->pd);
			res->pd = NULL;
		}
		if (res->ib_ctx)
		{
			ibv_close_device(res->ib_ctx);
			res->ib_ctx = NULL;
		}
		if (dev_list)
		{
			ibv_free_device_list(dev_list);
			dev_list = NULL;
		}
		if (res->sock >= 0)
		{
			if (close(res->sock))
				fprintf(stderr, "failed to close socket\n");
			res->sock = -1;
		}
	}
	return rc;
}
/******************************************************************************
* Function: modify_qp_to_init
*
* Input
* qp QP to transition
*
* Output
* none
*
* Returns
* 0 on success, ibv_modify_qp failure code on failure
*
* Description
* Transition a QP from the RESET to INIT state
******************************************************************************/
static int modify_qp_to_init(struct ibv_qp *qp)
{
	struct ibv_qp_attr attr;
	int flags;
	int rc;
	memset(&attr, 0, sizeof(attr));
	attr.qp_state = IBV_QPS_INIT;
	attr.port_num = config.ib_port;
	attr.pkey_index = 0;
	attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_ATOMIC;
	//attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE;
	flags = IBV_QP_STATE | IBV_QP_PKEY_INDEX | IBV_QP_PORT | IBV_QP_ACCESS_FLAGS;
	rc = ibv_modify_qp(qp, &attr, flags);
	if (rc)
		fprintf(stderr, "failed to modify QP state to INIT\n");
	return rc;
}
/******************************************************************************
* Function: modify_qp_to_rtr
*
* Input
* qp QP to transition
* remote_qpn remote QP number
* dlid destination LID
* dgid destination GID (mandatory for RoCEE)
*
* Output
* none
*
* Returns
* 0 on success, ibv_modify_qp failure code on failure
*
* Description
* Transition a QP from the INIT to RTR state, using the specified QP number
******************************************************************************/
static int modify_qp_to_rtr(struct ibv_qp *qp, uint32_t remote_qpn, uint16_t dlid, uint8_t *dgid)
{
	struct ibv_qp_attr attr;
	int flags;
	int rc;
	memset(&attr, 0, sizeof(attr));
	attr.qp_state = IBV_QPS_RTR;
	attr.path_mtu = IBV_MTU_256;
	attr.dest_qp_num = remote_qpn;
	attr.rq_psn = 0;
	attr.max_dest_rd_atomic = 1;
	attr.min_rnr_timer = 0x12;
	attr.ah_attr.is_global = 0;
	attr.ah_attr.dlid = dlid;
	attr.ah_attr.sl = 0;
	attr.ah_attr.src_path_bits = 0;
	attr.ah_attr.port_num = config.ib_port;
	if (config.gid_idx >= 0)
	{
		attr.ah_attr.is_global = 1;
		attr.ah_attr.port_num = 1;
		memcpy(&attr.ah_attr.grh.dgid, dgid, 16);
		attr.ah_attr.grh.flow_label = 0;
		attr.ah_attr.grh.hop_limit = 1;
		attr.ah_attr.grh.sgid_index = config.gid_idx;
		attr.ah_attr.grh.traffic_class = 0;
	}
	flags = IBV_QP_STATE | IBV_QP_AV | IBV_QP_PATH_MTU | IBV_QP_DEST_QPN |
			IBV_QP_RQ_PSN | IBV_QP_MAX_DEST_RD_ATOMIC | IBV_QP_MIN_RNR_TIMER;
	rc = ibv_modify_qp(qp, &attr, flags);
	if (rc)
		fprintf(stderr, "failed to modify QP state to RTR\n");
	return rc;
}
/******************************************************************************
* Function: modify_qp_to_rts
*
* Input
* qp QP to transition
*
* Output
* none
*
* Returns
* 0 on success, ibv_modify_qp failure code on failure
*
* Description
* Transition a QP from the RTR to RTS state
******************************************************************************/
static int modify_qp_to_rts(struct ibv_qp *qp)
{
	struct ibv_qp_attr attr;
	int flags;
	int rc;
	memset(&attr, 0, sizeof(attr));
	attr.qp_state = IBV_QPS_RTS;
	attr.timeout = 0x12;
	attr.retry_cnt = 6;
	attr.rnr_retry = 0;
	attr.sq_psn = 0;
	attr.max_rd_atomic = 1;
	flags = IBV_QP_STATE | IBV_QP_TIMEOUT | IBV_QP_RETRY_CNT |
			IBV_QP_RNR_RETRY | IBV_QP_SQ_PSN | IBV_QP_MAX_QP_RD_ATOMIC;
	rc = ibv_modify_qp(qp, &attr, flags);
	if (rc)
		fprintf(stderr, "failed to modify QP state to RTS\n");
	return rc;
}
/******************************************************************************
* Function: connect_qp
*
* Input
* res pointer to resources structure
*
* Output
* none
*
* Returns
* 0 on success, error code on failure
*
* Description
* Connect the QP. Transition the server side to RTR, sender side to RTS
******************************************************************************/
static int connect_qp(struct resources *res)
{	
	struct cm_con_data_t local_con_data;
	struct cm_con_data_t remote_con_data;
	struct cm_con_data_t tmp_con_data;
	int rc = 0;
	char temp_char;
	union ibv_gid my_gid;
	if (config.gid_idx >= 0)
	{
		rc = ibv_query_gid(res->ib_ctx, config.ib_port, config.gid_idx, &my_gid);
		if (rc)
		{
			fprintf(stderr, "could not get gid for port %d, index %d\n", config.ib_port, config.gid_idx);
			return rc;
		}
	}
	else
		memset(&my_gid, 0, sizeof my_gid);
	/* exchange using TCP sockets info required to connect QPs */
	local_con_data.addr = htonll((uintptr_t)res->buf);
	local_con_data.rkey = htonl(res->mr->rkey);
	local_con_data.qp_num = htonl(res->qp->qp_num);
	local_con_data.lid = htons(res->port_attr.lid);
	memcpy(local_con_data.gid, &my_gid, 16);
	fprintf(stdout, "\nLocal LID = 0x%x\n", res->port_attr.lid);
	if (sock_sync_data(res->sock, sizeof(struct cm_con_data_t), (char *)&local_con_data, (char *)&tmp_con_data) < 0)
	{
		fprintf(stderr, "failed to exchange connection data between sides\n");
		rc = 1;
		goto connect_qp_exit;
	}
	remote_con_data.addr = ntohll(tmp_con_data.addr);
	remote_con_data.rkey = ntohl(tmp_con_data.rkey);
	remote_con_data.qp_num = ntohl(tmp_con_data.qp_num);
	remote_con_data.lid = ntohs(tmp_con_data.lid);
	memcpy(remote_con_data.gid, tmp_con_data.gid, 16);
	/* save the remote side attributes, we will need it for the post SR */
	res->remote_props = remote_con_data;
	fprintf(stdout, "Remote address = 0x%" PRIx64 "\n", remote_con_data.addr);
	fprintf(stdout, "Remote rkey = 0x%x\n", remote_con_data.rkey);
	fprintf(stdout, "Remote QP number = 0x%x\n", remote_con_data.qp_num);
	fprintf(stdout, "Remote LID = 0x%x\n", remote_con_data.lid);
	if (config.gid_idx >= 0)
	{
		uint8_t *p = remote_con_data.gid;
		fprintf(stdout, "Remote GID =%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x:%02x\n ",p[0],
				  p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11], p[12], p[13], p[14], p[15]);
	}
	/* modify the QP to init */
	rc = modify_qp_to_init(res->qp);
	if (rc)
	{
		fprintf(stderr, "change QP state to INIT failed\n");
		goto connect_qp_exit;
	}
	/* let the client post RR to be prepared for incoming messages */
	// if (config.server_name)
	// {
	// 	rc = post_receive(res);
	// 	if (rc)
	// 	{
	// 		fprintf(stderr, "failed to post RR\n");
	// 		goto connect_qp_exit;
	// 	}
	// }
	rc = modify_qp_to_rtr(res->qp, remote_con_data.qp_num, remote_con_data.lid, remote_con_data.gid);
	if (rc)
	{
		fprintf(stderr, "failed to modify QP state to RTR\n");
		goto connect_qp_exit;
	}
	/* modify the QP to RTS */
	rc = modify_qp_to_rts(res->qp);
	if (rc)
	{
		fprintf(stderr, "failed to modify QP state to RTS\n");
		goto connect_qp_exit;
	}
	fprintf(stdout, "QP state was change to RTS\n");
	// rc = modify_qp_to_rtr(res->qp, remote_con_data.qp_num, remote_con_data.lid, remote_con_data.gid);
	// if (rc)
	// {
	// 	fprintf(stderr, "failed to modify QP state to RTR\n");
	// 	goto connect_qp_exit;
	// }
	/* sync to make sure that both sides are in states that they can connect to prevent packet loose */
	if (sock_sync_data(res->sock, 1, "Q", &temp_char)) /* just send a dummy char back and forth */
	{
		fprintf(stderr, "sync error after QPs are were moved to RTS\n");
		rc = 1;
	}
connect_qp_exit:
	return rc;
}
/******************************************************************************
* Function: resources_destroy
*
* Input
* res pointer to resources structure
*
* Output
* none
*
* Returns
* 0 on success, 1 on failure
*
* Description
* Cleanup and deallocate all resources used
******************************************************************************/
static int resources_destroy(struct resources *res)
{
	int rc = 0;
	if (res->qp)
		if (ibv_destroy_qp(res->qp))
		{
			fprintf(stderr, "failed to destroy QP\n");
			rc = 1;
		}
	if (res->mr)
		if (ibv_dereg_mr(res->mr))
		{
			fprintf(stderr, "failed to deregister MR\n");
			rc = 1;
		}
	if (res->buf)
		free(res->buf);
	if (res->cq)
		if (ibv_destroy_cq(res->cq))
		{
			fprintf(stderr, "failed to destroy CQ\n");
			rc = 1;
		}
	if (res->pd)
		if (ibv_dealloc_pd(res->pd))
		{
			fprintf(stderr, "failed to deallocate PD\n");
			rc = 1;
		}
	if (res->ib_ctx)
		if (ibv_close_device(res->ib_ctx))
		{
			fprintf(stderr, "failed to close device context\n");
			rc = 1;
		}
	if (res->sock >= 0)
		if (close(res->sock))
		{
			fprintf(stderr, "failed to close socket\n");
			rc = 1;
		}
	return rc;
}

void vranlc(struct resources *res, int n, double *x_seed, double a, uint64_t y){
	int i;
	double x,t1,t2,t3,t4,a1,a2,x1,x2,z;
	//printf("n = %d", n);
	/*
	 * ---------------------------------------------------------------------
	 * break A into two parts such that A = 2^23 * A1 + A2.
	 * ---------------------------------------------------------------------
	 */
	t1 = r23 * a;
	a1 = (int)t1;
	a2 = a - t23 * a1;
	x = *x_seed;

	/*
	 * ---------------------------------------------------------------------
	 * generate N results. this loop is not vectorizable.
	 * ---------------------------------------------------------------------
	 */
	//printf("n = %d\n", n);
	for(i=0; i<n; i++){
		/*
		 * ---------------------------------------------------------------------
		 * break X into two parts such that X = 2^23 * X1 + X2, compute
		 * Z = A1 * X2 + A2 * X1  (mod 2^23), and then
		 * X = 2^23 * Z + A2 * X2  (mod 2^46).
		 * ---------------------------------------------------------------------
		 */
		t1 = r23 * x;
		x1 = (int)t1;
		x2 = x - t23 * x1;
		t1 = a1 * x2 + a2 * x1;
		t2 = (int)(r23 * t1);
		z = t1 - t23 * t2;
		t3 = t23 * z + a2 * x2;
		t4 = (int)(r46 * t3);
		x = t3 - t46 * t4;
		//y[i] = r46 * x;
		res->buf[i] = r46 * x; 
		//printf("res buf[i] = %f\n", res->buf[i]);
		
	}
	mywrite(res, y, res->buf, sizeof(double) * n);
	if (poll_completion(res)) {
		fprintf(stderr, "poll completion failed\n");
	}
	//exit(0);
	*x_seed = x;
}

double randlc(double *x, double a){    
	double t1,t2,t3,t4,a1,a2,x1,x2,z;

	/*
	 * ---------------------------------------------------------------------
	 * break A into two parts such that A = 2^23 * A1 + A2.
	 * ---------------------------------------------------------------------
	 */
	t1 = r23 * a;
	a1 = (int)t1;
	a2 = a - t23 * a1;

	/*
	 * ---------------------------------------------------------------------
	 * break X into two parts such that X = 2^23 * X1 + X2, compute
	 * Z = A1 * X2 + A2 * X1  (mod 2^23), and then
	 * X = 2^23 * Z + A2 * X2  (mod 2^46).
	 * ---------------------------------------------------------------------
	 */
	t1 = r23 * (*x);
	x1 = (int)t1;
	x2 = (*x) - t23 * x1;
	t1 = a1 * x2 + a2 * x1;
	t2 = (int)(r23 * t1);
	z = t1 - t23 * t2;
	t3 = t23 * z + a2 * x2;
	t4 = (int)(r46 * t3);
	(*x) = t3 - t46 * t4;

	return (r46 * (*x));
}
