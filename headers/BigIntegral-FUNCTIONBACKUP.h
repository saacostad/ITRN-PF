
double calculate_f_NN_imaginary_integral(double q, double *params){

	double newParams[5] = {q, params[0], params[1], params[2], params[3]};

	
	// Creo una función lambda que tome una sola variable de integración)
	auto integrand = [&newParams](double b){
		return f_NN_imaginary_integrand(b, newParams);
	};


	double sum = integrand(a) + integrand(bMax);
	double x;

	#pragma omp parallel for reduction(+:sum)
	for (size_t i = 1; i < N; i++){
		
		x = a + (i * h);

		if (i % 2 == 0){
			sum += 2*integrand(x);
		} else {
			sum += 4*integrand(x);
		}

	}

	return sum * (h / 3.0);
}



double calculate_f_NN_real_integral(double q, double *params){

	double newParams[5] = {q, params[0], params[1], params[2], params[3]};

	
	// Creo una función lambda que tome una sola variable de integración)
	auto integrand = [&newParams](double b){
		return f_NN_real_integrand(b, newParams);
	};



	// INTEGRACIÓN POR SIMPSON PUES MARDITO GAUSS CONVERGE A 0

	// double sum = integrand(a) + integrand(bMax);
	double sum = f_NN_real_integrand(a, newParams) + f_NN_real_integrand(bMax, newParams);

	#pragma omp parallel for reduction(+:sum)
	for (size_t i = 1; i < N; i++){
		
		double x = a + (i * h);

		if (i % 2 == 0){
			// sum += 2*integrand(x);
			sum += 2* f_NN_real_integrand(x, newParams);
		} else {
			//  sum += 4*integrand(x);
			sum += 4*f_NN_real_integrand(x, newParams);
		}

	}

	return sum * (h / 3.0);
}


