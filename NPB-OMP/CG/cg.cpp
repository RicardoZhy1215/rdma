#include "omp.h"
#include "../common/npb-CPP.hpp"
#include "npbparams.hpp"
//#include "rdma.cpp"

#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <endian.h>
#include <byteswap.h>
#include <getopt.h>
// 或者
#include <malloc.h>
#include <atomic>
#include <sys/time.h>
#include <arpa/inet.h>
#include <infiniband/verbs.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include<algorithm>
#include <math.h>
#define MAX_POLL_CQ_TIMEOUT 2000

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
	double *buf;						   /* memory buffer pointer, used for RDMA and send ops */
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
int mywrite1(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size);
int myread1(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size);


const size_t block_num = 1000;
const size_t block_size = 2000;
const size_t cpupool_mem_size = 6020 * 1024 / 10;
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
           


/*
 * ---------------------------------------------------------------------
 * note: please observe that in the routine conj_grad three 
 * implementations of the sparse matrix-vector multiply have
 * been supplied. the default matrix-vector multiply is not
 * loop unrolled. the alternate implementations are unrolled
 * to a depth of 2 and unrolled to a depth of 8. please
 * experiment with these to find the fastest for your particular
 * architecture. if reporting timing results, any of these three may
 * be used without penalty.
 * ---------------------------------------------------------------------
 * class specific parameters: 
 * it appears here for reference only.
 * these are their values, however, this info is imported in the npbparams.h
 * include file, which is written by the sys/setparams.c program.
 * ---------------------------------------------------------------------
 */
#define NZ (NA*(NONZER+1)*(NONZER+1))
#define NAZ (NA*(NONZER+1))
#define T_INIT 0
#define T_BENCH 1
#define T_CONJ_GRAD 2
#define T_LAST 3

/* global variables */
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
static int colidx[NZ];
static int rowstr[NA+1];
static int iv[NA];
static int arow[NA];
static int acol[NAZ];
static double aelt[NAZ];
//static double a[NZ];
static double x[NA+2];
static double z[NA+2];
static double p[NA+2];
static double q[NA+2];
static double r[NA+2];
#else
static int (*colidx)=(int*)malloc(sizeof(int)*(NZ));
static int (*rowstr)=(int*)malloc(sizeof(int)*(NA+1));
static int (*iv)=(int*)malloc(sizeof(int)*(NA));
static int (*arow)=(int*)malloc(sizeof(int)*(NA));
static int (*acol)=(int*)malloc(sizeof(int)*(NAZ));
static double (*aelt)=(double*)malloc(sizeof(double)*(NAZ));

uint64_t a = 0;
static double (*x)=(double*)malloc(sizeof(double)*(NA+2));
static double (*z)=(double*)malloc(sizeof(double)*(NA+2));
static double (*p)=(double*)malloc(sizeof(double)*(NA+2));
static double (*q)=(double*)malloc(sizeof(double)*(NA+2));
static double (*r)=(double*)malloc(sizeof(double)*(NA+2));
#endif
static int naa;
static int nzz;
static int firstrow;
static int lastrow;
static int firstcol;
static int lastcol;
static double amult;
static double tran;
static boolean timeron;

/* function prototypes */
static void conj_grad(struct resources *res,int colidx[],
		int rowstr[],
		double x[],
		double z[],
		uint64_t a,
		double p[],
		double q[],
		double r[],
		double* rnorm);
static int icnvrt(double x,
		int ipwr2);
static void makea(struct resources *res,int n,
		int nz,
		uint64_t a,
		int colidx[],
		int rowstr[],
		int firstrow,
		int lastrow,
		int firstcol,
		int lastcol,
		int arow[],
		int acol[][NONZER+1],
		double aelt[][NONZER+1],
		int iv[]);
static void sparse(struct resources *res, uint64_t a,
		int colidx[],
		int rowstr[],
		int n,
		int nz,
		int nozer,
		int arow[],
		int acol[][NONZER+1],
		double aelt[][NONZER+1],
		int firstrow,
		int lastrow,
		int nzloc[],
		double rcond,
		double shift);
static void sprnvc(int n,
		int nz,
		int nn1,
		double v[],
		int iv[]);
static void vecset(int n,
		double v[],
		int iv[],
		int* nzv,
		int i,
		double val);

/* cg */
int main(int argc, char **argv){
#if defined(DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION)
	printf(" DO_NOT_ALLOCATE_ARRAYS_WITH_DYNAMIC_MEMORY_AND_AS_SINGLE_DIMENSION mode on\n");
#endif
	
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
	u_int64_t slices_cnt = ceil(sizeof(double) * NZ / cpupool_mem_size) + 1;
	printf("slices is %d\n", slices_cnt);
	if (sizeof(double) * NZ < cpupool_mem_size || sizeof(double) * NZ == cpupool_mem_size) {
		a = mymalloc(&res, sizeof(double) * NZ, 0);
		if (poll_completion(&res)) {
			fprintf(stderr, "poll completion failed\n");
		}
	} else {
		for (int i = 0; i < slices_cnt; i++) {
			if (i == slices_cnt - 1) {
				addr = mymalloc(&res, sizeof(double) * NZ - (slices_cnt - 1) * cpupool_mem_size, i * cpupool_mem_size);
				if (i == 0) {
					a = addr;
				}
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
			} else {
				addr = mymalloc(&res, cpupool_mem_size, i * cpupool_mem_size); 
				if (i == 0) {
					a = addr;
				}
				if (poll_completion(&res)) {
					fprintf(stderr, "poll completion failed\n");
				}
			}
		}
	}

	int	i, j, k, it;
	double zeta;
	double rnorm;
	double norm_temp1, norm_temp2;
	double t, mflops, tmax;
	char class_npb;
	boolean verified;
	double zeta_verify_value, epsilon, err;

	char *t_names[T_LAST];

	for(i=0; i<T_LAST; i++){
		timer_clear(i);
	}

	FILE* fp;
	if((fp = fopen("timer.flag", "r")) != NULL){
		timeron = TRUE;
		t_names[T_INIT] = (char*)"init";
		t_names[T_BENCH] = (char*)"benchmk";
		t_names[T_CONJ_GRAD] = (char*)"conjgd";
		fclose(fp);
	}else{
		timeron = FALSE;
	}

	timer_start(T_INIT);

	firstrow = 0;
	lastrow  = NA-1;
	firstcol = 0;
	lastcol  = NA-1;

	if(NA == 1400 && NONZER == 7 && NITER == 15 && SHIFT == 10.0){
		class_npb = 'S';
		zeta_verify_value = 8.5971775078648;
	}else if(NA == 7000 && NONZER == 8 && NITER == 15 && SHIFT == 12.0){
		class_npb = 'W';
		zeta_verify_value = 10.362595087124;
	}else if(NA == 14000 && NONZER == 11 && NITER == 15 && SHIFT == 20.0){
		class_npb = 'A';
		zeta_verify_value = 17.130235054029;
	}else if(NA == 75000 && NONZER == 13 && NITER == 75 && SHIFT == 60.0){
		class_npb = 'B';
		zeta_verify_value = 22.712745482631;
	}else if(NA == 150000 && NONZER == 15 && NITER == 75 && SHIFT == 110.0){
		class_npb = 'C';
		zeta_verify_value = 28.973605592845;
	}else if(NA == 1500000 && NONZER == 21 && NITER == 100 && SHIFT == 500.0){
		class_npb = 'D';
		zeta_verify_value = 52.514532105794;
	}else if(NA == 9000000 && NONZER == 26 && NITER == 100 && SHIFT == 1500.0){
		class_npb = 'E';
		zeta_verify_value = 77.522164599383;
	}else{
		class_npb = 'U';
	}

	printf("\n\n NAS Parallel Benchmarks 4.1 Parallel C++ version with OpenMP - CG Benchmark\n\n");
	printf(" Size: %11d\n", NA);
	printf(" Iterations: %5d\n", NITER);

	naa = NA;
	nzz = NZ;

	/* initialize random number generator */
	tran    = 314159265.0;
	amult   = 1220703125.0;
	zeta    = randlc( &tran, amult );
		//13049872 12431376
	makea(&res, naa, 
			nzz, 
			a, 
			colidx, 
			rowstr, 
			firstrow, 
			lastrow, 
			firstcol, 
			lastcol, 
			arow, 
			(int(*)[NONZER+1])(void*)acol, 
			(double(*)[NONZER+1])(void*)aelt,
			iv);



	#pragma omp parallel private(it,i,j,k)	
	{
		#pragma omp for nowait

		for(j = 0; j < lastrow - firstrow + 1; j++){
			for(k = rowstr[j]; k < rowstr[j+1]; k++){
				colidx[k] = colidx[k] - firstcol;
			}
		}
		/* set starting vector to (1, 1, .... 1) */
		#pragma omp for nowait
		for(i = 0; i < NA+1; i++){
			x[i] = 1.0;
		}
		#pragma omp for nowait
		for(j = 0; j<lastcol-firstcol+1; j++){
			q[j] = 0.0;
			z[j] = 0.0;
			r[j] = 0.0;
			p[j] = 0.0;
		}
		
		#pragma omp single
			zeta = 0.0;

		/*
		 * -------------------------------------------------------------------
		 * ---->
		 * do one iteration untimed to init all code and data page tables
		 * ----> (then reinit, start timing, to niter its)
		 * -------------------------------------------------------------------*/
	
		for(it = 1; it <= 1; it++){
			/* the call to the conjugate gradient routine */

			conj_grad(&res, colidx, rowstr, x, z, a, p, q, r, &rnorm);
			#pragma omp single
			{
				norm_temp1 = 0.0;
				norm_temp2 = 0.0;
			}
			
			/*
			 * --------------------------------------------------------------------
			 * zeta = shift + 1/(x.z)
			 * so, first: (x.z)
			 * also, find norm of z
			 * so, first: (z.z)
			 * --------------------------------------------------------------------
			 */
			#pragma omp for reduction(+:norm_temp1,norm_temp2)
			for(j = 0; j < lastcol - firstcol + 1; j++){
				norm_temp1 += x[j] * z[j];
				norm_temp2 += + z[j] * z[j];
			}

			#pragma omp single
				norm_temp2 = 1.0 / sqrt(norm_temp2);

			/* normalize z to obtain x */
			#pragma omp for
			for(j = 0; j < lastcol - firstcol + 1; j++){     
				x[j] = norm_temp2 * z[j];
			}

		} /* end of do one iteration untimed */
	// FILE *file_a;
    // file_a = fopen("a1.txt", "w"); 
    // if (file_a == NULL) {
    //     perror("Error opening file");
    // }
	// int cnt = ceil(NZ * sizeof(double) / cpupool_mem_size) + 1;
	// printf("cnt is %d\n", cnt);
	// for (int i = 0; i < cnt; i++) {
	// 	if (i == cnt - 1) {
	// 		myread(&res, a + i * cpupool_mem_size, res.buf, NZ * sizeof(double) - i * cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	} else {
	// 		myread(&res, a + i * cpupool_mem_size, res.buf, cpupool_mem_size);
	// 		if (poll_completion(&res)) {
	// 			fprintf(stderr, "poll completion failed\n");
	// 		}
	// 	}
	// 	for (int j = 0; j < std::min(NZ * sizeof(double) - i * cpupool_mem_size, cpupool_mem_size) / sizeof(double); j++) {
	// 		fprintf(file_a, "%f \n", res.buf[j]);
	// 	}
	// }
	// exit(0);
		// double profile_result = double(profile_end - profile_start) / CLOCKS_PER_SEC;
		// printf( "poll time is %f seconds\n", duration);
		// printf( "RDMA time is %f seconds\n", rdma_result);
		// printf( "profile result is %f seconds\n", profile_result);
		
		/* set starting vector to (1, 1, .... 1) */	
		#pragma omp for
		for(i = 0; i < NA+1; i++){
			x[i] = 1.0;
		}

		#pragma omp single
			zeta = 0.0;

		#pragma omp master
		{
			timer_stop(T_INIT);

			printf(" Initialization time = %15.3f seconds\n", timer_read(T_INIT));
			
			timer_start(T_BENCH);
		}
		
		// printf("exit here!\n");
		// exit(0);
		/*
		 * --------------------------------------------------------------------
		 * ---->
		 * main iteration for inverse power method
		 * ---->
		 * --------------------------------------------------------------------
		 */
		double loop_start, loop_end;
		
		loop_start = omp_get_wtime();
		for(it = 1; it <= NITER; it++){
			/* the call to the conjugate gradient routine */
			#pragma omp master
			if(timeron){timer_start(T_CONJ_GRAD);}
			conj_grad(&res, colidx, rowstr, x, z, a, p, q, r, &rnorm);
			#pragma omp master
			if(timeron){timer_stop(T_CONJ_GRAD);}

			#pragma omp single
			{
				norm_temp1 = 0.0;
				norm_temp2 = 0.0;
			}

			/*
			 * --------------------------------------------------------------------
			 * zeta = shift + 1/(x.z)
			 * so, first: (x.z)
			 * also, find norm of z
			 * so, first: (z.z)
			 * --------------------------------------------------------------------
			 */
			#pragma omp for reduction(+:norm_temp1,norm_temp2)
			for(j = 0; j < lastcol - firstcol + 1; j++){
				norm_temp1 += x[j]*z[j];
				norm_temp2 += z[j]*z[j];
			}
			#pragma omp single
			{
				norm_temp2 = 1.0 / sqrt(norm_temp2);
				zeta = SHIFT + 1.0 / norm_temp1;
			}

			#pragma omp master
			{
				if(it==1){printf("\n   iteration           ||r||                 zeta\n");}
				printf("    %5d       %20.14e%20.13e\n", it, rnorm, zeta);
			}
			/* normalize z to obtain x */
			#pragma omp for 
			for(j = 0; j < lastcol - firstcol + 1; j++){
				x[j] = norm_temp2 * z[j];
			}

			

		} /* end of main iter inv pow meth */

		loop_end = omp_get_wtime();
		//printf("loop timer is: %f\n", loop_end - loop_start);

	} /* end parallel */
	timer_stop(T_BENCH);


	/*
	 * --------------------------------------------------------------------
	 * end of timed section
	 * --------------------------------------------------------------------
	 */

	t = timer_read(T_BENCH);

	if (sock_sync_data(res.sock, 1, "R", &temp_char)) /* just send a dummy char back and forth */
	{
		fprintf(stderr, "sync error before RDMA ops\n");
		rc = 1;
	}

	printf(" Benchmark completed\n");

	epsilon = 1.0e-10;
	if(class_npb != 'U'){
		err = fabs(zeta - zeta_verify_value) / zeta_verify_value;
		if(err <= epsilon){
			verified = TRUE;
			printf(" VERIFICATION SUCCESSFUL\n");
			printf(" Zeta is    %20.13e\n", zeta);
			printf(" Error is   %20.13e\n", err);
		}else{
			verified = FALSE;
			printf(" VERIFICATION FAILED\n");
			printf(" Zeta                %20.13e\n", zeta);
			printf(" The correct zeta is %20.13e\n", zeta_verify_value);
		}
	}else{
		verified = FALSE;
		printf(" Problem size unknown\n");
		printf(" NO VERIFICATION PERFORMED\n");
	}
	if(t != 0.0){
		mflops = (double)(2.0*NITER*NA)
			* (3.0+(double)(NONZER*(NONZER+1))
					+ 25.0
					* (5.0+(double)(NONZER*(NONZER+1)))+3.0)
			/ t / 1000000.0;
	}else{
		mflops = 0.0;
	}
	setenv("OMP_NUM_THREADS","1",0);
	c_print_results((char*)"CG",
			class_npb,
			NA,
			0,
			0,
			NITER,
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
		if(tmax == 0.0){tmax = 1.0;}
		printf("  SECTION   Time (secs)\n");
		for(i = 0; i < T_LAST; i++){
			t = timer_read(i);
			if(i == T_INIT){
				printf("  %8s:%9.3f\n", t_names[i], t);
			}else{
				printf("  %8s:%9.3f  (%6.2f%%)\n", t_names[i], t, t*100.0/tmax);
				if(i == T_CONJ_GRAD){
					t = tmax - t;
					printf("    --> %8s:%9.3f  (%6.2f%%)\n", "rest", t, t*100.0/tmax);
				}
			}
		}
	}
	return 0;
}

/*
 * ---------------------------------------------------------------------
 * floating point arrays here are named as in NPB1 spec discussion of 
 * CG algorithm
 * ---------------------------------------------------------------------
 */
static void conj_grad(struct resources *res, int colidx[],
		int rowstr[],
		double x[],
		double z[],
		uint64_t a,
		double p[],
		double q[],
		double r[],
		double* rnorm){
	int j, k;
	int cgit, cgitmax;
	double alpha, beta, suml;
	static double d, sum, rho, rho0;

	cgitmax = 25;
	#pragma omp single nowait
	{

		rho = 0.0;
		sum = 0.0;
	}
	/* initialize the CG algorithm */
	#pragma omp for
	for(j = 0; j < naa+1; j++){
		q[j] = 0.0;
		z[j] = 0.0;
		r[j] = x[j];
		p[j] = r[j];
	}
 
	/*
	 * --------------------------------------------------------------------
	 * rho = r.r
	 * now, obtain the norm of r: First, sum squares of r elements locally...
	 * --------------------------------------------------------------------
	 */
	#pragma omp for reduction(+:rho)
	for(j = 0; j < lastcol - firstcol + 1; j++){
		rho += r[j]*r[j];
	}

	

	/* the conj grad iteration loop */
	for(cgit = 1; cgit <= cgitmax; cgit++){
		/*
		 * ---------------------------------------------------------------------
		 * q = A.p
		 * the partition submatrix-vector multiply: use workspace w
		 * ---------------------------------------------------------------------
		 * 
		 * note: this version of the multiply is actually (slightly: maybe %5) 
		 * faster on the sp2 on 16 nodes than is the unrolled-by-2 version 
		 * below. on the Cray t3d, the reverse is TRUE, i.e., the 
		 * unrolled-by-two version is some 10% faster.  
		 * the unrolled-by-8 version below is significantly faster
		 * on the Cray t3d - overall speed of code is 1.5 times faster.
		 */

		#pragma omp single nowait
		{
			d = 0.0;
			/*
			 * --------------------------------------------------------------------
			 * save a temporary of rho
			 * --------------------------------------------------------------------
			 */
			rho0 = rho;
			rho = 0.0;
		}
		#pragma omp for nowait
		for(j = 0; j < lastrow - firstrow + 1; j++){
			suml = 0.0;
			int tid = omp_get_thread_num();
			int x = 0;
			int n = 0;
			struct timespec req, rem;
			req.tv_sec = 0;       
    		req.tv_nsec = 1;     
			if(sizeof(double) * (rowstr[j+1] - rowstr[j]) <= block_size ) {
				myread(res, a + rowstr[j] * sizeof(double), res->buf + tid * block_size, sizeof(double) * (rowstr[j+1] - rowstr[j]));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
			
			} 
			// else if (sizeof(double) * (rowstr[j+1] - rowstr[j]) > block_size && sizeof(double) * (rowstr[j+1] - rowstr[j]) < cpupool_mem_size) {
			// 	myread(res, a + rowstr[j] * sizeof(double), res->buf, sizeof(double) * (rowstr[j+1] - rowstr[j]));
			// 	if (poll_completion(res)) {
			// 		fprintf(stderr, "poll completion failed\n");
			// 	}
			// 	printf("hello\n");
			// }
			//  else {// degrade to single thread
			// 	int loop_cnt = ceil(sizeof(double) * (rowstr[j+1] - rowstr[j]) / (block_size)) + 1;
			// 	for (int i = 0; i < loop_cnt; i++) {
			// 		if (i == loop_cnt - 1) {
			// 			myread(res, a + k * sizeof(double) + i * block_size, res->buf, sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]) - i * block_size);
			// 			//Todo
			// 		} else {
			// 			myread(res, a + k * sizeof(double) + i * block_size, res->buf, block_size);
			// 			//
			// 		}
				
			// 	}
			// }

			// if (nanosleep(&req, &rem) < 0) {
        	// 	perror("Nanosleep failed");
    		// }

			for(k = rowstr[j]; k < rowstr[j+1]; k++){
				//suml += /*a[k]*/ res->buf[/*tid*/ /* block_size +*/ tid * block_size + i] * p[colidx[k]];
				suml += res->buf[tid * block_size  + x] * p[colidx[k]];
				x++;
			}
			q[j] = suml;
		}
		/*timer*/
		/*
		 * --------------------------------------------------------------------
		 * obtain p.q
		 * --------------------------------------------------------------------
		 */

		#pragma omp for reduction(+:d)
		for (j = 0; j < lastcol - firstcol + 1; j++) {
			d += p[j]*q[j];
		}

		/*
		 * --------------------------------------------------------------------
		 * obtain alpha = rho / (p.q)
		 * -------------------------------------------------------------------
		 */
		alpha = rho0 / d;
			
		/*
		 * ---------------------------------------------------------------------
		 * obtain z = z + alpha*p
		 * and    r = r - alpha*q
		 * ---------------------------------------------------------------------
		 */

		#pragma omp for reduction(+:rho)
		for(j = 0; j < lastcol - firstcol + 1; j++){
			z[j] += alpha*p[j];
			r[j] -= alpha*q[j];

			/*
			 * ---------------------------------------------------------------------
			 * rho = r.r
			 * now, obtain the norm of r: first, sum squares of r elements locally...
			 * ---------------------------------------------------------------------
			 */
			rho += r[j]*r[j];
		}

		/*
		 * ---------------------------------------------------------------------
		 * obtain beta
		 * ---------------------------------------------------------------------
		 */	
		beta = rho / rho0;

		/*
		 * ---------------------------------------------------------------------
		 * p = r + beta*p
		 * ---------------------------------------------------------------------
		 */
		#pragma omp for
		for(j = 0; j < lastcol - firstcol + 1; j++){
			p[j] = r[j] + beta*p[j];
		}
	} /* end of do cgit=1, cgitmax */




	/*
	 * ---------------------------------------------------------------------
	 * compute residual norm explicitly: ||r|| = ||x - A.z||
	 * first, form A.z
	 * the partition submatrix-vector multiply
	 * ---------------------------------------------------------------------
	 */
	#pragma omp for nowait
	for(j = 0; j < lastrow - firstrow + 1; j++){
		suml = 0.0;
		int tid = omp_get_thread_num();
		int m = 0;
		if(sizeof(double) * (rowstr[j+1] - rowstr[j]) <= block_size) {
			myread(res, a + rowstr[j] * sizeof(double), res->buf + tid * block_size, sizeof(double) * (rowstr[j+1] - rowstr[j]));
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
		} 
			// else if (sizeof(double) * (rowstr[j+1] - rowstr[j]) > block_size && sizeof(double) * (rowstr[j+1] - rowstr[j]) < cpupool_mem_size) {
			// 	myread(res, a + rowstr[j] * sizeof(double), res->buf, sizeof(double) * (rowstr[j+1] - rowstr[j]));
			// 	if (poll_completion(res)) {
			// 		fprintf(stderr, "poll completion failed\n");
			// 	}
			// }
			//  else {// degrade to single thread
			// 	int loop_cnt = ceil(sizeof(double) * (rowstr[j+1] - rowstr[j]) / block_size) + 1;
			// 	for (int i = 0; i < loop_cnt; i++) {
			// 		if (i == loop_cnt - 1) {
			// 			myread(res, a + k * sizeof(double) + i * block_size, res->buf, sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]) - i * block_size);
			// 			//Todo
			// 		} else {
			// 			myread(res, a + k * sizeof(double) + i * block_size, res->buf, block_size);
			// 			//
			// 		}
				
			// 	}
			// }
		// myread(res, a[0] + rowstr[j] * sizeof(double), res->buf + tid * block_size, sizeof(double) * (rowstr[j+1] - rowstr[j]));
		// if (poll_completion(res)) {
		// 	fprintf(stderr, "poll completion failed\n"); 
		// }

		for(k = rowstr[j]; k < rowstr[j+1] ; k++){
			suml += /*a[k]*/ res->buf[tid * block_size + m/*tid*/] * z[colidx[k]];
			m++;
		}
		r[j] = suml;
	}


	/*
	 * ---------------------------------------------------------------------
	 * at this point, r contains A.z
	 * ---------------------------------------------------------------------
	 */
	#pragma omp for reduction(+:sum)
	for(j = 0; j < lastcol-firstcol+1; j++){
		suml   = x[j] - r[j];
		sum += suml*suml;
	}
	#pragma omp single
		*rnorm = sqrt(sum);
}

/*
 * ---------------------------------------------------------------------
 * scale a double precision number x in (0,1) by a power of 2 and chop it
 * ---------------------------------------------------------------------
 */
static int icnvrt(double x, int ipwr2){
	return (int)(ipwr2 * x);
}

/*
 * ---------------------------------------------------------------------
 * generate the test problem for benchmark 6
 * makea generates a sparse matrix with a
 * prescribed sparsity distribution
 *
 * parameter    type        usage
 *
 * input
 *
 * n            i           number of cols/rows of matrix
 * nz           i           nonzeros as declared array size
 * rcond        r*8         condition number
 * shift        r*8         main diagonal shift
 *
 * output
 *
 * a            r*8         array for nonzeros
 * colidx       i           col indices
 * rowstr       i           row pointers
 *
 * workspace
 *
 * iv, arow, acol i
 * aelt           r*8
 * ---------------------------------------------------------------------
 */
static void makea(struct resources *res, int n,
		int nz,
		uint64_t a,
		int colidx[],
		int rowstr[],
		int firstrow,
		int lastrow,
		int firstcol,
		int lastcol,
		int arow[],
		int acol[][NONZER+1],
		double aelt[][NONZER+1],
		int iv[]){
	int iouter, ivelt, nzv, nn1;
	int ivc[NONZER+1];

	double vc[NONZER+1];
	/*
	 * --------------------------------------------------------------------
	 * nonzer is approximately  (int(sqrt(nnza /n)));
	 * --------------------------------------------------------------------
	 * nn1 is the smallest power of two not less than n
	 * --------------------------------------------------------------------
	 */
	nn1 = 1;
	do{
		nn1 = 2 * nn1;
	}while(nn1 < n);

	/*
	 * -------------------------------------------------------------------
	 * generate nonzero positions and save for the use in sparse
	 * -------------------------------------------------------------------
	 */


	for(iouter = 0; iouter < n; iouter++){
		nzv = NONZER;
		sprnvc(n, nzv, nn1, vc, ivc);
		vecset(n, vc, ivc, &nzv, iouter+1, 0.5);
		arow[iouter] = nzv;
		for(ivelt = 0; ivelt < nzv; ivelt++){
			acol[iouter][ivelt] = ivc[ivelt] - 1;
			aelt[iouter][ivelt] = vc[ivelt];
		}
	}
	/*
	 * ---------------------------------------------------------------------
	 * ... make the sparse matrix from list of elements with duplicates
	 * (iv is used as  workspace)
	 * ---------------------------------------------------------------------
	 */
	sparse(res, a,
			colidx,
			rowstr,
			n,
			nz,
			NONZER,
			arow,
			acol,
			aelt,
			firstrow,
			lastrow,
			iv,
			RCOND,
			SHIFT);
}

/*
 * ---------------------------------------------------------------------
 * rows range from firstrow to lastrow
 * the rowstr pointers are defined for nrows = lastrow-firstrow+1 values
 * ---------------------------------------------------------------------
 */
static void sparse(struct resources *res, uint64_t a,
		int colidx[],
		int rowstr[],
		int n,
		int nz,
		int nozer,
		int arow[],
		int acol[][NONZER+1],
		double aelt[][NONZER+1],
		int firstrow,
		int lastrow,
		int nzloc[],
		double rcond,
		double shift){	
	int nrows;
	char temp_char;

	/*
	 * ---------------------------------------------------
	 * generate a sparse matrix from a list of
	 * [col, row, element] tri
	 * ---------------------------------------------------
	 */
	int i, j, j1, j2, nza, k, kk, nzrow, jcol;
	double size, scale, ratio, va;
	double sum;
	boolean goto_40;
	/*
	 * --------------------------------------------------------------------
	 * how many rows of result
	 * --------------------------------------------------------------------
	 */
	nrows = lastrow - firstrow + 1;

	/*
	 * --------------------------------------------------------------------
	 * ...count the number of triples in each row
	 * --------------------------------------------------------------------
	 */
	for(j = 0; j < nrows+1; j++){
		rowstr[j] = 0;
	}
	for(i = 0; i < n; i++){
		for(nza = 0; nza < arow[i] ; nza++){
			j = acol[i][nza] + 1;
			rowstr[j] = rowstr[j] + arow[i];
		}
	}
	rowstr[0] = 0;
	for(j = 1; j < nrows+1; j++){
		rowstr[j] = rowstr[j] + rowstr[j-1];
	}
	nza = rowstr[nrows] - 1;

	/*
	 * ---------------------------------------------------------------------
	 * ... rowstr(j) now is the location of the first nonzero
	 * of row j of a
	 * ---------------------------------------------------------------------
	 */
	if(nza > nz){
		printf("Space for matrix elements exceeded in sparse\n");
		printf("nza, nzmax = %d, %d\n", nza, nz);
		exit(EXIT_FAILURE);
	}

	/*
	 * ---------------------------------------------------------------------
	 * ... preload data pages
	 * ---------------------------------------------------------------------
	 */

	/**timer1************************************************************************/
	//clock_t start1 = clock();

	for(j = 0; j < nrows; j++){
		if (sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]) <= cpupool_mem_size) {
			memset(res->buf, 0 , sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]));
			mywrite(res, a + k * sizeof(double), res->buf, sizeof(double) * (rowstr[j+1] - rowstr[j]));
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
		} 
		// else {
		// 	int loop_cnt = ceil(sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]) / cpupool_mem_size) + 1;
		// 	memset(res->buf, 0 , cpupool_mem_size);
		// 	for (int i = 0; i < loop_cnt; i++) {
		// 		if (i == loop_cnt - 1) {
		// 			mywrite(res, a + k * sizeof(double) + i * cpupool_mem_size, res->buf, sizeof(res->buf[0]) * (rowstr[j+1] - rowstr[j]) - i * cpupool_mem_size);
		// 			if (poll_completion(res)) {
		// 				fprintf(stderr, "poll completion failed\n");
		// 			}
		// 		} else {
		// 			mywrite(res, a + k * sizeof(double) + i * cpupool_mem_size, res->buf, cpupool_mem_size);
		// 			if (poll_completion(res)) {
		// 				fprintf(stderr, "poll completion failed\n");
		// 			}
		// 		}
				
		// 	}
		// }


		for(k = rowstr[j]; k < rowstr[j+1]; k++){
			// res->buf[0] = 0.0;
			// mywrite(res, a[0] + k * sizeof(double), res->buf + 0, sizeof(res->buf[0]));
			// if (poll_completion(res)) {
			// 	fprintf(stderr, "poll completion failed\n");
			// }
			//a[k] = 0.0;
			colidx[k] = -1;
		}
		nzloc[j] = 0;
	}

	
	size = 1.0;
	ratio = pow(rcond, (1.0 / (double)(n)));

	clock_t start2 = clock();
	for(i = 0; i < n; i++){
		for(nza = 0; nza < arow[i] ; nza++){
			j = acol[i][nza];
			scale = size * aelt[i][nza];
			for(nzrow = 0; nzrow < arow[i]; nzrow++){
				jcol = acol[i][nzrow];
				va = aelt[i][nzrow] * scale;
		
				if(jcol == j && j == i){
					va = va + rcond - shift;
				}

				goto_40 = FALSE;
				/**timer2************************************************************************/
				// int index = 0;
				// k = rowstr[j];
				// myread(res, a + k * sizeof(double), res->buf, sizeof(double) * ((rowstr[j+1] - 2) - k + 1));
				// if (poll_completion(res)) {
				// 	fprintf(stderr, "poll completion failed\n");
				// }
				for(k = rowstr[j]; k < rowstr[j+1]; k++){
					if(colidx[k] > jcol){
						if (sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1) <= cpupool_mem_size) {
                            //clock_t test_start = clock();
							myread(res, a + k * sizeof(double), res->buf, sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
							
							mywrite(res, a + (k + 1) * sizeof(double), res->buf, sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1));
							if (poll_completion(res)) {
								fprintf(stderr, "poll completion failed\n");
							}
                            // clock_t test_end = clock();
							// sum +=  (double) (test_end - test_start) / CLOCKS_PER_SEC;
							// memcpy(res->buf[!index], res->buf[index], sizeof(double) * ((rowstr[j+1] - 2) - k + 1));
							
							// mywrite(res, a + (k + 1) * sizeof(double), res->buf[!index], sizeof(double) * ((rowstr[j+1] - 2) - k + 1));
							// if (poll_completion(res)) {
							// 	fprintf(stderr, "poll completion failed\n");
							// }
							// myread(res, a + (k + 1) * sizeof(double), res->buf[index], sizeof(double) * ((rowstr[j+1] - 2) - k + 1));
							// if (poll_completion(res)) {
							// 	fprintf(stderr, "poll completion failed\n");
							// }
						} 
						// else {
						// 	int loop_cnt = ceil(sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1)) / cpupool_mem_size + 1;
						// 	for (int i = 0; i < loop_cnt; i++) {
						// 		if (i == loop_cnt - 1) {
						// 			myread(res, a + k * sizeof(double) + i * cpupool_mem_size, res->buf, sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1) - cpupool_mem_size);
						// 			if (poll_completion(res)) {
						// 				fprintf(stderr, "poll completion failed\n");
						// 			}
						// 			mywrite(res, a + (k + 1) * sizeof(double) + i * cpupool_mem_size, res->buf, sizeof(res->buf[0]) * ((rowstr[j+1] - 2) - k + 1) - cpupool_mem_size);
						// 			if (poll_completion(res)) {
						// 				fprintf(stderr, "poll completion failed\n");
						// 			}
						// 		} else {
						// 			myread(res, a + k * sizeof(double) + i * cpupool_mem_size, res->buf, cpupool_mem_size);
						// 			if (poll_completion(res)) {
						// 				fprintf(stderr, "poll completion failed\n");
						// 			}
						// 			mywrite(res, a + (k + 1) * sizeof(double) + i * cpupool_mem_size, res->buf, cpupool_mem_size);
						// 			if (poll_completion(res)) {
						// 				fprintf(stderr, "poll completion failed\n");
						// 			}
						// 		}
						// 	}	
						// }
						for(kk = rowstr[j+1]-2; kk >= k; kk--){
							if(colidx[kk] > -1){
								colidx[kk+1] = colidx[kk];
							}
						}
						colidx[k] = jcol;
						res->buf[0] = 0.0;
						mywrite(res, a + k * sizeof(double), res->buf + 0, sizeof(res->buf[0]));
						if (poll_completion(res)) {
							fprintf(stderr, "poll completion failed\n");
						}
						//a[k]  = 0.0;
						goto_40 = TRUE;
						break;
					}else if(colidx[k] == -1){
						colidx[k] = jcol;
						goto_40 = TRUE;
						break;
					}else if(colidx[k] == jcol){
						/*
						 * --------------------------------------------------------------
						 * ... mark the duplicated entry
						 * -------------------------------------------------------------
						 */
						nzloc[j] = nzloc[j] + 1;
						goto_40 = TRUE;
						break;
					}
				}
				
				

				if(goto_40 == FALSE){
					printf("internal error in sparse: i=%d\n", i);
					exit(EXIT_FAILURE);
				}
				myread(res, a + k * sizeof(double), res->buf + 0, sizeof(res->buf[0]));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				res->buf[0] += va;
				mywrite(res, a + k * sizeof(double), res->buf + 0, sizeof(res->buf[0]));
				if (poll_completion(res)) {
					fprintf(stderr, "poll completion failed\n");
				}
				//a[k] = a[k] + va;
			}
		}
		size = size * ratio;
	}
	// clock_t end2 = clock();
	// printf("timer for sparse is: %lf\n", (double) (end2 - start2) / CLOCKS_PER_SEC);
    // printf("timer for rdma is: %lf\n", sum);

	/*
	 * ---------------------------------------------------------------------
	 * ... remove empty entries and generate final results
	 * ---------------------------------------------------------------------
	 */
	for(j = 1; j < nrows; j++){
		nzloc[j] = nzloc[j] + nzloc[j-1];
	}

	/**timer3************************************************************************/
	clock_t start3 = clock();
	for(j = 0; j < nrows; j++){
		if(j > 0){
			j1 = rowstr[j] - nzloc[j-1];
		}else{
			j1 = 0;
		}
		j2 = rowstr[j+1] - nzloc[j];
		nza = rowstr[j];
		for(k = j1; k < j2; k++){
			myread(res, a + nza * sizeof(double), res->buf + 0, sizeof(res->buf[0]));
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}

			mywrite(res, a + k * sizeof(double),  res->buf + 0, sizeof(res->buf[0]));
			if (poll_completion(res)) {
				fprintf(stderr, "poll completion failed\n");
			}
			//a[k] = a[nza];
			colidx[k] = colidx[nza];
			nza = nza + 1;
		}
	}
	clock_t end3 = clock();

	for(j = 1; j < nrows+1; j++){
		rowstr[j] = rowstr[j] - nzloc[j-1];
	}
	nza = rowstr[nrows] - 1;
	
}

/*
 * ---------------------------------------------------------------------
 * generate a sparse n-vector (v, iv)
 * having nzv nonzeros
 *
 * mark(i) is set to 1 if position i is nonzero.
 * mark is all zero on entry and is reset to all zero before exit
 * this corrects a performance bug found by John G. Lewis, caused by
 * reinitialization of mark on every one of the n calls to sprnvc
 * ---------------------------------------------------------------------
 */
static void sprnvc(int n, int nz, int nn1, double v[], int iv[]){
	int nzv, ii, i;
	double vecelt, vecloc;

	nzv = 0;

	while(nzv < nz){
		vecelt = randlc(&tran, amult);

		/*
		 * --------------------------------------------------------------------
		 * generate an integer between 1 and n in a portable manner
		 * --------------------------------------------------------------------
		 */
		vecloc = randlc(&tran, amult);
		i = icnvrt(vecloc, nn1) + 1;
		if(i>n){continue;}

		/*
		 * --------------------------------------------------------------------
		 * was this integer generated already?
		 * --------------------------------------------------------------------
		 */
		boolean was_gen = FALSE;
		for(ii = 0; ii < nzv; ii++){
			if(iv[ii] == i){
				was_gen = TRUE;
				break;
			}
		}
		if(was_gen){continue;}
		v[nzv] = vecelt;
		iv[nzv] = i;
		nzv = nzv + 1;
	}
}

/*
 * --------------------------------------------------------------------
 * set ith element of sparse vector (v, iv) with
 * nzv nonzeros to val
 * --------------------------------------------------------------------
 */
static void vecset(int n, double v[], int iv[], int* nzv, int i, double val){
	int k;
	boolean set;

	set = FALSE;
	for(k = 0; k < *nzv; k++){
		if(iv[k] == i){
			v[k] = val;
			set  = TRUE;
		}
	}
	if(set == FALSE){
		v[*nzv]  = val;
		iv[*nzv] = i;
		*nzv     = *nzv + 1;
	}
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
// /******************************************************************************
// End of socket operations
// ******************************************************************************/
// /* poll_completion */
// /******************************************************************************
// * Function: poll_completion
// *
// * Input
// * res pointer to resources structure
// *
// * Output
// * none
// *
// * Returns
// * 0 on success, 1 on failure
// *
// * Description
// * Poll the completion queue for a single event. This function will continue to
// * poll the queue until MAX_POLL_CQ_TIMEOUT milliseconds have passed.
// *
// ******************************************************************************/
static int poll_completion(struct resources *res)
{
	struct ibv_wc wc[1];
	// unsigned long start_time_msec;
	// unsigned long cur_time_msec;
	// struct timeval cur_time;
	int poll_result;
	int rc = 0;

	do
	{
		poll_result = ibv_poll_cq(res->cq, 500, wc);

	} while (poll_result == 0/*((cur_time_msec - start_time_msec) < MAX_POLL_CQ_TIMEOUT)*/);
	if (poll_result < 0) {
		/* poll CQ failed */
		fprintf(stderr, "poll CQ failed\n");
		rc = 1;
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
	clock_t poll2_start, poll2_end;
	double poll2_res = 0;
	poll2_start = clock();
	do
	{
	poll_result = ibv_poll_cq(res->cq, 1, &wc);
	} while ((poll_result == 0));
	poll2_end = clock();
	poll2_res += (double)(poll2_end - poll2_start)/ CLOCKS_PER_SEC;
	//printf( "poll2_Res time is %f seconds\n",  poll2_res);
	if (poll_result < 0)
	{
		/* poll CQ failed */
		fprintf(stderr, "poll CQ failed\n");
		rc = 1;
	}
	// else if (poll_result == 0)
	// { /* the CQ is empty */
	// 	fprintf(stderr, "completion wasn't found in the CQ after timeout\n");
	// 	rc = 1;
	// }
	// else
	// {
	// 	/* CQE found */
	// 	//fprintf(stdout, "completion was found in CQ with status 0x%x\n", wc.status);
	// 	/* check the completion status (here we don't care about the completion opcode */
	// 	if (wc.status != IBV_WC_SUCCESS)
	// 	{
	// 		fprintf(stderr, "got bad completion with status: 0x%x, vendor syndrome: 0x%x\n", wc.status, wc.vendor_err);
	// 		exit(0);
	// 		rc = 1;
	// 	}
	// }
	return rc;
}


// // static int post_malloc(struct resources *res, size_t size)
// // {
// // 	struct ibv_send_wr sr;
// // 	struct ibv_sge sge;
// // 	struct ibv_send_wr *bad_wr = NULL;
// // 	int rc;
// // 	/* prepare the scatter/gather entry */
// // 	memset(&sge, 0, sizeof(sge));
// // 	sge.addr = (uintptr_t)res->buf;
// // 	sge.length = block_size;
// // 	sge.lkey = res->mr->lkey;
// // 	/* prepare the send work request */
// // 	memset(&sr, 0, sizeof(sr));
// // 	sr.next = NULL;
// // 	sr.wr_id = 0;
// // 	sr.sg_list = &sge;
// // 	sr.num_sge = 1;
// // 	sr.opcode = IBV_WR_SEND;
// // 	sr.send_flags = IBV_SEND_SIGNALED;
// // 	if (opcode == IBV_WR_ATOMIC_CMP_AND_SWP) {
// //         // 原子比较并交换操作
// //         sr.wr.atomic.remote_addr = res->remote_props.addr;
// //         sr.wr.atomic.rkey = res->remote_props.rkey;
// //         sr.wr.atomic.compare_add = res->atomic_compare;
// //         sr.wr.atomic.swap = res->atomic_swap;
// //     } else if (opcode == IBV_WR_ATOMIC_FETCH_AND_ADD) {
// //         // 原子取并加操作
// //         sr.wr.atomic.remote_addr = res->remote_props.addr;
// //         sr.wr.atomic.rkey = res->remote_props.rkey;
// //         sr.wr.atomic.compare_add = res->atomic_add;
// // 	} else if (opcode != IBV_WR_SEND) {
// // 		sr.wr.rdma.remote_addr = res->remote_props.addr;
// //         sr.wr.rdma.rkey = res->remote_props.rkey;
// // 	}

// // 	/* there is a Receive Request in the responder side, so we won't get any into RNR flow */
// // 	rc = ibv_post_send(res->qp, &sr, &bad_wr);
// // 	if (rc) {
// // 		fprintf(stderr, "failed1 to post SR\n");
// // 	}
// //     else {
// // 		printf(stdout, "Send Request was posted\n");
// // 	}
// // 		return rc;
// // 	}
	



// /******************************************************************************
// * Function: mymalloc
// ******************************************************************************/

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

int myread1(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size) {
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
int mywrite1(struct resources *res, uint64_t remote_mem_addr, double* local_buf, size_t size) {
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
// /******************************************************************************
// * Function: post_receive
// *
// * Input
// * res pointer to resources structure
// *
// * Output
// * none
// *
// * Returns
// * 0 on success, error code on failure
// *
// * Description
// *
// ******************************************************************************/
// static int post_receive(struct resources *res)
// {
// 	struct ibv_recv_wr rr;
// 	struct ibv_sge sge;
// 	struct ibv_recv_wr *bad_wr;
// 	int rc;
// 	/* prepare the scatter/gather entry */
// 	memset(&sge, 0, sizeof(sge));
// 	sge.addr = (uintptr_t)res->buf;
// 	sge.length = 1024;
// 	sge.lkey = res->mr->lkey;
// 	/* prepare the receive work request */
// 	memset(&rr, 0, sizeof(rr));
// 	rr.next = NULL;
// 	rr.wr_id = 0;
// 	rr.sg_list = &sge;
// 	rr.num_sge = 1;
// 	/* post the Receive Request to the RQ */
// 	rc = ibv_post_recv(res->qp, &rr, &bad_wr);
// 	if (rc)
// 		fprintf(stderr, "failed to post RR\n");
// 	else
// 		fprintf(stdout, "Receive Request was posted\n");
// 	return rc;
// }
// /******************************************************************************
// * Function: resources_init
// *
// * Input
// * res pointer to resources structure
// *
// * Output
// * res is initialized
// *
// * Returns
// * none
// *
// * Description
// * res is initialized to default values
// ******************************************************************************/
static void resources_init(struct resources *res)
{
	memset(res, 0, sizeof *res);
	res->sock = -1;
}
// /******************************************************************************
// * Function: resources_create
// *
// * Input
// * res pointer to resources structure to be filled in
// *
// * Output
// * res filled in with resources
// *
// * Returns
// * 0 on success, 1 on failure
// *
// * Description
// *
// * This function creates and allocates all necessary system resources. These
// * are stored in res.
// *****************************************************************************/
static int resources_create(struct resources *res)
{
	struct ibv_device **dev_list = NULL;
	struct ibv_qp_init_attr qp_init_attr;
	struct ibv_device *ib_dev = NULL;
	//struct ibv_device_attr device_attr;
	size_t size;
	int i;
	int mr_flags = 0;
	int cq_size = 500;
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
	//memset(res->buf1, 0.0, size);
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
	res->mr= ibv_reg_mr(res->pd, res->buf, size, mr_flags);
	//res->mr[1] = ibv_reg_mr(res->pd, res->buf[1], size, mr_flags);
	if (!res->mr)
	{
		fprintf(stderr, "ibv_reg_mr failed with mr_flags=0x%x\n", mr_flags);
		rc = 1;
		goto resources_create_exit;
	}
	// if (!res->mr1)
	// {
	// 	fprintf(stderr, "ibv_reg_mr failed with mr_flags=0x%x\n", mr_flags);
	// 	rc = 1;
	// 	goto resources_create_exit;
	// }
	fprintf(stdout, "MR was registered with addr=%p, lkey=0x%x, rkey=0x%x, flags=0x%x\n",
			res->buf, res->mr->lkey, res->mr->rkey, mr_flags);
	/* create the Queue Pair */
	memset(&qp_init_attr, 0, sizeof(qp_init_attr));
	qp_init_attr.qp_type = IBV_QPT_RC;
	qp_init_attr.sq_sig_all = 1;
	qp_init_attr.send_cq = res->cq;
	qp_init_attr.recv_cq = res->cq;
	qp_init_attr.cap.max_send_wr = 1000;
	qp_init_attr.cap.max_recv_wr = 100;
	qp_init_attr.cap.max_send_sge = 10;
	qp_init_attr.cap.max_recv_sge = 10;
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
// /******************************************************************************
// * Function: modify_qp_to_rtr
// *
// * Input
// * qp QP to transition
// * remote_qpn remote QP number
// * dlid destination LID
// * dgid destination GID (mandatory for RoCEE)
// *
// * Output
// * none
// *
// * Returns
// * 0 on success, ibv_modify_qp failure code on failure
// *
// * Description
// * Transition a QP from the INIT to RTR state, using the specified QP number
// ******************************************************************************/
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
// /******************************************************************************
// * Function: connect_qp
// *
// * Input
// * res pointer to resources structure
// *
// * Output
// * none
// *
// * Returns
// * 0 on success, error code on failure
// *
// * Description
// * Connect the QP. Transition the server side to RTR, sender side to RTS
// ******************************************************************************/
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
// /******************************************************************************
// * Function: resources_destroy
// *
// * Input
// * res pointer to resources structure
// *
// * Output
// * none
// *
// * Returns
// * 0 on success, 1 on failure
// *
// * Description
// * Cleanup and deallocate all resources used
// ******************************************************************************/
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
