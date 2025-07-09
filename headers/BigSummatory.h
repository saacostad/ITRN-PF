#pragma once

#include <cmath>
#include <iostream>
#include <complex>
#include <omp.h>					// Para computación en paralelo
#include <chrono>

#include "BigIntegral.h"				// Donde se encuentran las definiciones de las funciones del Eikonal
#include "AuxiliarFunction.h"				// Donde se encuentran las definiciones de Coulomb
#include <boost/math/special_functions/legendre.hpp>	// Para los pol. legendre


extern double K;
double inK = 1.0 / K;

using complex = std::complex<double>;
using std::exp;




int lMax = 300;

//==================================================================
//	CORRIMIENTO DE FASE TOTAL
//==================================================================
complex Sel(double b, double *params){
    return exp( complex(0.0, 1.0) * calculate_XI_integral(b, params) );
}


//==================================================================
//	DEFINICIÓN DE LA bl
//==================================================================
double bl(double l){
    return inK * ( eta + std::sqrt( (eta * eta) + ((l + 0.5)*(l + 0.5)) ) ) ;
}



//==================================================================
//	DEFINICIÓN DE LA AMPLITUD DE SCATTERING TOTAL Fel
//==================================================================
//

complex FelSummand(int l, double theta, double *params){
    return (2.0*double(l) + 1.0) * exp( 2.0 * complex(0.0, 1.0) * CoulombPhaseShift(double(l))) * (1.0 - Sel( bl(double(l)), params )) * boost::math::legendre_p(l, std::cos(theta));
}

complex Fel(double theta, double *params){

    // return CoulombScatteringAmplitude(theta);
    complex Fc = CoulombScatteringAmplitude(theta) + CoulombScatteringAmplitude(pi - theta);
    complex constant(0.0, 1.0 / (K));

    
    double realSum = 0.0;
    double imagSum = 0.0;

    auto start2 = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for reduction(+:realSum, imagSum)
    for (int l = 0; l < lMax; l += 2){

	complex val = FelSummand(l, theta, params); 
	
	realSum += val.real();
	imagSum += val.imag();
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = end2 - start2;
//    std::cout << "Elapsed time: " << elapsed2.count() << " seconds\n"
//    std::cout << "Sumatory: " << realSum << " + " << imagSum << "i" << std::endl;
//    std::cout << "Constant: " << constant << std::endl;
//   std::cout << "Fc: " << Fc << std::endl;
    
    return complex(realSum, imagSum) * constant + Fc;
}





//==================================================================
//	PRUEBA DEL HEADER
//==================================================================
void runSumHeader(void){
    double params[4] = {-281.5, 0.0, 0.99, 1};
    std::cout << Fel(0.1, params) << std::endl;
    

    //std::cout << CoulombScatteringAmplitude(0.1) << std::endl;
    //std::cout << FelSummand(1, 0.1, params) << std::endl;

    return;
}
