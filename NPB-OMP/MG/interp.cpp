int index = 0;
		myread(res, z + (0 * mm2 * mm1) * sizeof(double), res->buf + index * mid, sizeof(double) * 2 * mm2 * mm1);
		myread(res, u + (2 * 0 * n2 * n1) * sizeof(double), res->buf + index * mid + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
		#pragma omp for
		for(i3 = 0; i3 < mm3-1; i3++){
			int tid = 0; 
			int num_element = min(2 * mm2 * mm1 + 2 * n1 * n2, static_cast<int>(block_size / sizeof(double)));
				for (int i = 0; i < 1; i++) {
					if (poll_completion_num(res, 2)) {
						fprintf(stderr, "poll completion failed\n");
					}
						
				for(i2 = 0; i2 < mm2-1; i2++){
				for(i1 = 0; i1 < mm1; i1++){
					z1[i1] = res->buf[index * mid + tid * block_size / sizeof(double) + (i2 + 1) * mm1 + i1] + res->buf[index * mid + tid * block_size / sizeof(double) + i2 * mm1 + i1];
					z2[i1] = res->buf[index * mid + tid * block_size / sizeof(double) + mm2 * mm1 + i2 * mm1 + i1] + res->buf[index * mid + tid * block_size / sizeof(double) + i2 * mm1 + i1];
					z3[i1] = res->buf[index * mid + tid * block_size / sizeof(double) + mm2 * mm1 + (i2 + 1) * mm1 + i1] + res->buf[index * mid + tid * block_size / sizeof(double) + mm2 * mm1 + i2 * mm1 + i1] + z1[i1];
					// z1[i1] = z[i3][i2+1][i1] + z[i3][i2][i1];
					// z2[i1] = z[i3+1][i2][i1] + z[i3][i2][i1];
					// z3[i1] = z[i3+1][i2+1][i1] + z[i3+1][i2][i1] + z1[i1];
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1] += res->buf[index * mid + tid * block_size / sizeof(double) + i2 * mm1 + i1];
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + 2 * i2 * n1 + 2 * i1 + 1] +=  0.5 * (res->buf[index * mid + tid * block_size / sizeof(double) + i2 * mm1 + i1 + 1] + res->buf[index * mid + tid * block_size / sizeof(double) + i2 * mm1 + i1]);
					// u[2*i3][2*i2][2*i1] = u[2*i3][2*i2][2*i1] + z[i3][i2][i1];
					// u[2*i3][2*i2][2*i1+1] = u[2*i3][2*i2][2*i1+1] + 0.5 * (z[i3][i2][i1+1]+z[i3][i2][i1]);
				}
				
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1] += 0.5 * z1[i1];
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.25 * (z1[i1] + z1[i1 + 1]);
					// u[2*i3][2*i2+1][2*i1] = u[2*i3][2*i2+1][2*i1] + 0.5 * z1[i1];
					// u[2*i3][2*i2+1][2*i1+1] = u[2*i3][2*i2+1][2*i1+1] + 0.25 * ( z1[i1] + z1[i1+1] );
					
				}

				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1] += 0.5 * z2[i1];
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + 2 * i2 * n1 + 2 * i1 + 1] += 0.25 * (z2[i1] + z2[i1 + 1]);

					// u[2*i3+1][2*i2][2*i1] = u[2*i3+1][2*i2][2*i1] + 0.5 * z2[i1];
					// u[2*i3+1][2*i2][2*i1+1] = u[2*i3+1][2*i2][2*i1+1]+ 0.25*( z2[i1] + z2[i1+1] );
					
				}
				for(i1 = 0; i1 < mm1-1; i1++){
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1] += 0.25 * z3[i1];
					res->buf[index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1 + 1] += 0.125 * (z3[i1] + z3[i1 + 1]);
					// u[2*i3+1][2*i2+1][2*i1] = u[2*i3+1][2*i2+1][2*i1] + 0.25* z3[i1];
					// u[2*i3+1][2*i2+1][2*i1+1] = u[2*i3+1][2*i2+1][2*i1+1] + 0.125*( z3[i1] + z3[i1+1] );
					//printf("%f %f\n", res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1], res->buf[tid * block_size / sizeof(double) + 2 * mm2 * mm1 + n1 * n2 + (2 * i2 + 1) * n1 + 2 * i1 + 1]);
				}
			}
			if (i3 != 0) {
				if (poll_completion_num(res, 2)) {
					fprintf(stderr, "poll completion failed\n");
				}
			}	
			if (i3 + 1 >= mm3-1) {
				goto skip_code;
			}  
			myread(res, z + ((i3 + 1) * mm2 * mm1 + i * num_element) * sizeof(double), res->buf + (!index) * mid + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
			myread(res, u + (2 * (i3 + 1) * n2 * n1 + i * num_element) * sizeof(double), res->buf + (!index) * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);
		skip_code:
			mywrite(res, z + (i3 * mm2 * mm1 + i * num_element) * sizeof(double), res->buf + index * mid + tid * block_size / sizeof(double), sizeof(double) * 2 * mm2 * mm1);
			mywrite(res, u + (2 * i3 * n2 * n1 + i * num_element) * sizeof(double), res->buf + index * mid + tid * block_size / sizeof(double) + 2 * mm2 * mm1, sizeof(double) * 2 * n2 * n1);

			if (i3 == mm3 - 2) {
				if (poll_completion_num(res, 2)) {
					fprintf(stderr, "poll completion failed\n");
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
	 	index = !index;
	}
	clock_t end = clock();
	time = (double) (end - start) / CLOCKS_PER_SEC;
	time_sum += time;
