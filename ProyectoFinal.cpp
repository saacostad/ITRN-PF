#include <iostream>					// La confiable IOSTREAM
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>					// Para matemáticas en general



//-------------------------------------------------------
//		CONSTANTES DEL SISTEMA 
//-------------------------------------------------------

#include "headers/CONSTANTS.h"				// Para las constantes del sistema
double E_NN_lab = 85.0;

// double M = mN * A;							// Masa de cada núcleo
// double E_lab = E_NN_lab * A;						// Energía total del lab
// double E_cm = std::sqrt( std::pow( 2 * M , 2) + (2 * M  * E_lab) );	// Esto es directamente considerando la E. relativista
// double K = std::sqrt( std::pow(E_cm,2) - std::pow(2*M,2) ) / 2.0;	// Momento en el centro de masa
// double v_cm = K / std::sqrt( K*K + M*M );				// Velocidad en el centro de masa
// double E_tot_NN = mN + E_NN_lab;					// E total por nucleón
// double k = std::sqrt( E_tot_NN * E_tot_NN - mN * mN );			// Momento por nucleón
// double eta = (Z * Z * alpha * mu) / (K);			// Parámetro de Sommerfield
// double AB = A * A;
// double waveNumber = std::sqrt(2 * mu * E_cm) / hbar;


double E_lab = E_NN_lab * A;
double E_cm = E_lab * 0.5;
double K = std::sqrt(2.0 * mu * E_cm ) / hbar;
double v_cm = (K * hbar) / mu;

double k = std::sqrt(2 * mN * E_lab) / hbar;

double eta = (Z * Z * e2) / (hbar * v_cm);

double waveNumber = std::sqrt(2 * mu * E_cm) / hbar;
double AB = A * A;


//------------------------------------------------------	
//		CONSTANTES DE LAS FUNCIONES 
//------------------------------------------------------	


#include "headers/AuxiliarFunction.h"			// En donde están las funciones de Coulomb y el bl
double xi_NN_constant = -1.0 / (v_cm * hbar);


// INCLUSIÓN DE LAS LIBRERÍAS PROPIAS
#include "headers/BigIntegral.h"			// En donde se encuentran las integrales feas
#include "headers/BigSummatory.h"			// En este header se define el Sel y el Fel



double CrossSection(double theta, double *params){
	return std::norm( Fel(theta, params) / CoulombScatteringAmplitude(theta));
}




int main(void){
	omp_set_nested(1);
	// double params[4] = {-281.5, 0.0, 0.99, 1};
	double params[4] = {-67.5, 689.6, 0.72, 2.55}; // E_nn_lab = 85.0

	std::cerr << v_cm << std::endl;	
	for (double i = 1; i <= 14; i += 0.05){
		std::cerr << "Iteración " << i << std::endl;
		std::cout << double(i) << "\t" << CrossSection(double(i * pi / 180.0), params) << std::endl;
	}


	return 0;
}
