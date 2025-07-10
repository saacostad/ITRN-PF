#pragma once

#include <cmath>		// Como en este header voy a agregar las funciones de los factores de forma,
				// y para ello necesito la exponencial, entonces la agrego aquí

double Z = 6.0;			// No de protones
double A = 12.0;			// No másico   
//
//
double hbar = 197.327;		// hbar c [MeV fm]
double e2 = 1.44;		// carga eléctrica	
double mN = 931.5;		// Masa del nucleón

double mu = 5589.0;		// Masa reducida

double alpha = 1.0 / 137.036;


double formFactor(double q){
	return 1.494 * std::exp(- (q*q) * 0.741 ) - 0.494 * std::exp(- (q*q) * 0.417 );
}
