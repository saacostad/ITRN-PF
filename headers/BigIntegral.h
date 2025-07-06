#pragma once

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>	// Librerías para las funciones especiales
#include <boost/math/quadrature/gauss_kronrod.hpp>	// Para hacer la integral de f_NN
#include <boost/math/policies/policy.hpp>



#include <chrono>



double pi = 3.141592653589793238462643383279502884; 	// Para no importar una librería interea para el pi
extern double xi_NN_constant;				// Constantes para las funciones



double bMax = 8.0; 		// Tomamos este bMax como límite de integración 
int N = 50000;			// No intervales en el método de Simpson
double a = 1e-5;			// Limite inferior
double h = (bMax - a) / N;		// Cortes



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
void calculate_f_NN_integral_Gauss(double q, double* result, double *params){
	
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

	
	result[0] = real_result; result[1] = imaginary_result;
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
	for (int i = 1; i < 1000; i++){

		// double integral_val[2];
		// calculate_f_NN_integral_Gauss(i / 500, integral_val, params);

		calculate_f_NN_imaginary_integral_withGauss(i / 500, params);
		calculate_f_NN_imaginary_integral_withGauss(i / 500, params);
		if (i % 250 == 0){
			// std::cout << "Parte real: " << integral_val[0] << std::endl;
			//std::cout << "Parte imaginaria: " << integral_val[1] << std::endl;
		}
	}

	auto end2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - start2;




	std::cout << "Elapsed time using Gauss: " << elapsed2.count() << " seconds\n";




	return 0;
}


