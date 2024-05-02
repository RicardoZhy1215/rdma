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
	clock_t start = clock();
	double time = 0;
	int i3, i2, i1, d1, d2, d3, t1, t2, t3;

	double z1[M], z2[M], z3[M];
	if(timeron){
		#pragma omp master
			timer_start(T_INTERP);
	}

	if(n1 != 3 && n2 != 3 && n3 != 3){
		int index = 0;
		myread(res, z + (0 * mm2 * mm1) * sizeof(double), res->buf + index * mid, sizeof(double) * 2 * mm2 * mm1);
		myread(res, u + (2 * 0 * n2 * n1) * sizeof(double), res->buf + index * mid + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
		#pragma omp for
		for(i3 = 0; i3 < mm3-1; i3++){
			if (poll_completion_num(res, 2)) {
				fprintf(stderr, "poll completion failed\n");
			}		
			for(i2 = 0; i2 < mm2-1; i2++){

				if (i2==0 && i3 != 0) {
					if (poll_completion_num(res, 2)) {
						fprintf(stderr, "poll completion failed\n");
					}
				}

				if (i2 == 0 && i3 + 1 < mm3-1) {
					myread(res, z + ((i3 + 1) * mm2 * mm1 ) * sizeof(double), res->buf + (!index) * mid , sizeof(double) * 2 * mm2 * mm1);
					myread(res, u + (2 * (i3 + 1) * n2 * n1 ) * sizeof(double), res->buf + (!index) * mid + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
				}
				
				for(i1 = 0; i1 < mm1; i1++){
					z1[i1] = res->buf[index * mid + (i2 + 1) * mm1 + i1] + res->buf[index * mid + i2 * mm1 + i1];
					z2[i1] = res->buf[index * mid + mm2 * mm1 + i2 * mm1 + i1] + res->buf[index * mid + i2 * mm1 + i1];
					z3[i1] = res->buf[index * mid + mm2 * mm1 + (i2 + 1) * mm1 + i1] + res->buf[index * mid + mm2 * mm1 + i2 * mm1 + i1] + z1[i1];
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1] += res->buf[index * mid + i2 * mm1 + i1];
					res->buf[index * mid + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1 + 1] +=  0.5 * (res->buf[index * mid  + i2 * mm1 + i1 + 1] + res->buf[index * mid + i2 * mm1 + i1]);
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1] += 0.5 * z1[i1];
					res->buf[index * mid + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.25 * (z1[i1] + z1[i1 + 1]);
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1] += 0.5 * z2[i1];
					res->buf[index * mid + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1 + 1] += 0.25 * (z2[i1] + z2[i1 + 1]);
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1] += 0.25 * z3[i1];
					res->buf[index * mid + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.125 * (z3[i1] + z3[i1 + 1]);
				}
			}
	
			mywrite(res, z + (i3 * mm2 * mm1) * sizeof(double), res->buf + index * mid , sizeof(double) * 2 * mm2 * mm1);
			mywrite(res, u + (2 * i3 * n2 * n1) * sizeof(double), res->buf + index * mid + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
	 		index = !index;
		}
		if (poll_completion_num(res, 2)) {
			fprintf(stderr, "poll completion failed\n");
		}
		clock_t end = clock();
		time = (double) (end - start) / CLOCKS_PER_SEC;
		time_sum += time;
	}
