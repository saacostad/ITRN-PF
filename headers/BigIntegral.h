#pragma once

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>	// Librerías para las funciones especiales
#include <boost/math/quadrature/gauss_kronrod.hpp>	// Para hacer la integral de f_NN
#include <boost/math/policies/policy.hpp>
#include <omp.h>


#include <chrono>



double pi = 3.141592653589793238462643383279502884; 	// Para no importar una librería interea para el pi
extern double xi_NN_constant;				// Constantes para las funciones
std::complex<double> i(0.0, 1.0);			// Unidad imaginaria.
extern double AB;


				// bMax usual 6.0
double bMax = 5.5; 		// Tomamos este bMax como límite de integración en f_NN 

double qMax = 4.0;		// Tomamos este qMax como límite de integración en Xi
double a = 1e-6;
int N = 800;			// Intervalos en los que partir el grid en Simpson

//===============================================================================
//		POTENCIAL A AJUSTAR
//===============================================================================


// Función de la integral \int_{-\infty}^{\infty} V_NN(\sqrt(b^2 + z^2)) dz, 
// que gracias a dios es analítica
double V_NN_integral(double b, double V1, double V2, double a1, double a2){
	using namespace std;
	return (V1 * exp(-(a1 * a1) * (b*b) ) * (sqrt(pi) / a1))
	     + (V2 * exp(-(a2 * a2) * (b*b) ) * (sqrt(pi) / a2));
}


//=================================================================================
//	CORRIMIENTO DE FASE EIKONAL XI_{NN}
//=================================================================================


// Corrimiento de fase eikonal por nucleón
double Xi_NN(double b, double *params){
	return xi_NN_constant * V_NN_integral(b, params[0], params[1], params[2], params[3]); 	
}



//=================================================================================
//	INTEGRANDOS PARA LA AMPLITUD DE SCATERING EIKONAL F_NN
//=================================================================================


// Integrando de la parte real de la amplitud de scattering por nucleón
double f_NN_real_integrand(double b, double *params){
	return std::cyl_bessel_j(0, params[0] * b) * b * 
	       (1.0 - std::cos( Xi_NN(b, &params[1]) ) );
}

// Integrando de la parte imaginaria de la amplitud de scattering por nucleón
double f_NN_imaginary_integrand(double b, double *params){
	return -std::cyl_bessel_j(0, params[0] * b) * b * 
	       ( std::sin( Xi_NN(b, &params[1]) ) );
}




//=================================================================================
//	CÁLCULO DE LAS INTEGRALES PARA F_NN
//=================================================================================
// NOTA: las integrales en la librería no pueden correr en complejos, por lo tanto,
//me toca separarlas en sus componentes reales e imaginarios


// Cálculo de la integral f_NN para su componente real
std::complex<double> calculate_f_NN_integral_Gauss(double q, double *params){
	
	// Creo un nuevo arreglo de parámetros (para poder tener una única variable de integración)
	// que va a tener tanto los parámetros del potencial como q en su primer elemento
	double newParams[5] = {q, params[0], params[1], params[2], params[3]};

	
	// Creo una función lambda que tome una sola variable de integración
	auto real_integrand = [&newParams](double b){
		return f_NN_real_integrand(b, newParams);
		// Para la parte real
	};

	auto imaginary_integrand = [&newParams](double b){
		return f_NN_imaginary_integrand(b, newParams);
		// Para la parte imaginaria
	};



	// REALIZAMOS LAS INTEGRALES
	//_______________________________________________________________________

	double real_abs_error;
	double real_result = boost::math::quadrature::gauss_kronrod<double, 21, boost::math::policies::policy<> >::integrate(
	    real_integrand, 0.0, bMax,		// Función, bMin, bMax
	    12,					// Máximos subintervalos
	    1e-8,				// Tolerancia
	    &real_abs_error,
	    nullptr     
	);

	
	// Lo mismo pero para la parte imaginaria
	double imaginary_abs_error;
	double imaginary_result = boost::math::quadrature::gauss_kronrod<double, 21, boost::math::policies::policy<> >::integrate(
	    imaginary_integrand, 0.0, bMax,		// Función, bMin, bMax
	    12,					// Máximos subintervalos
	    1e-8,				// Tolerancia
	    &imaginary_abs_error,	
	    nullptr     
	);


	// Como f_NN es esta integral multiplicada por i*k (donde podemos ignorar la k pues desaparece en el \Xi), 
	// simplemente a nuestro resultado r + ib, lo multiplicamos por i para llegar a -b + imaginary_result
	// result[0] = -imaginary_result; result[1] = real_result;
	return std::complex<double>(real_result, imaginary_result) * i;
}




//==========================================================================
//		INTEGRAL PARA LA XI TOTAL
//==========================================================================


// Integrando de la Xi total
std::complex<double> XI_integrand(double q, double b, double *params){
	return q * std::pow( formFactor(q) , 2) * std::cyl_bessel_j(0, q * b) *
	       calculate_f_NN_integral_Gauss(q, params);
}


// Aquí, vamos a hallar Xi total para un b dado, haciendo la 
// integral por el método de Simpson
std::complex<double> calculate_XI_integral(double b, double *params){

	double h = (qMax - a) / N;		// Intervalos a dividir la integral
	
	double realSumOdd = 0.0; double imagSumOdd = 0.0;
	double realSumEven = 0.0; double imagSumEven = 0.0;
	

	#pragma omp parallel for reduction(+:realSumOdd, imagSumOdd)
	for (int i = 1; i < N; i +=2){
		std::complex<double> val = XI_integrand(a + (i*h), b, params);

		realSumOdd += val.real();
		imagSumOdd += val.imag();
	}


	#pragma omp parallel for reduction(+:realSumEven, imagSumEven)
	for (int i = 2; i < N; i +=2){
		std::complex<double> val = XI_integrand(a + (i*h), b, params);

		realSumEven += val.real();
		imagSumEven += val.imag();
	}

	
	std::complex<double> sumOdd(realSumOdd, imagSumOdd);
	std::complex<double> sumEven(realSumEven, imagSumEven);


	return AB * (h / 3.0) * (XI_integrand(a, b, params) + XI_integrand(qMax, b, params) + (4.0 * sumOdd) + (2.0 * sumEven));
}

	



//==========================================================================
//		PARA HACER PRUEBAS CON LAS FUNCIONES DE ESTE HEADER
//==========================================================================

double runHeader(void){
	double params[4] = {-281.5, 0.0, 0.99, 1};

//	for (int i = 1; i < 12; i++){
//		std::cout << std::endl;
//		std::cout << calculate_f_NN_imaginary_integral(0.4 * i, params) << std::endl;
//		std::cout << calculate_f_NN_real_integral(0.4 * i, params) << std::endl;
//		std::cout << std::endl;
//	}

	auto start2 = std::chrono::high_resolution_clock::now();

	
	#pragma omp parallel for 
	for (int i = 1; i < 125; i++){

		std::complex<double> val = calculate_XI_integral(0.1, params);
		std::cout << "qMax:" << qMax << " || N: " << N << std::endl;
		std::cout << val.real() << "+ " << val.imag() << "i" << std::endl;
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - start2;




	std::cout << "Elapsed time using Gauss: " << elapsed2.count() << " seconds\n";




	return 0;
}


