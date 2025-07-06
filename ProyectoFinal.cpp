#include <iostream>					// La confiable IOSTREAM
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>					// Para matemáticas en general

// INCLUSIÓN DE LAS LIBRERÍAS PROPIAS
#include "headers/CONSTANTS.h"				// Para las constantes del sistema
#include "headers/BigIntegral.h"



//-------------------------------------------------------
//		CONSTANTES DEL SISTEMA 
//-------------------------------------------------------

double E_NN_lab = 25.0;

double M = mN * A;							// Masa de cada núcleo
double E_lab = E_NN_lab * A;						// Energía total del lab
double E_cm = std::sqrt( std::pow( 2 * M , 2) + (2 * M  * E_lab) );	// Esto es directamente considerando la E. relativista
double K = std::sqrt( std::pow(E_cm,2) - std::pow(2*M,2) ) / 2.0;	// Momento en el centro de masa
double v_cm = K / std::sqrt( K*K + M*M );				// Velocidad en el centro de masa
double E_tot_NN = mN + E_NN_lab;					// E total por nucleón
double k = std::sqrt( E_tot_NN * E_tot_NN - mN * mN );			// Momento por nucleón
double eta = (Z * Z * (1.0 / 137.0) * mu) / (K);			// Parámetro de Sommerfield



//------------------------------------------------------	
//		CONSTANTES DE LAS FUNCIONES 
//------------------------------------------------------	


double xi_NN_constant = -1.0 / (v_cm);

int main(void){
	
	runHeader();
	return 0;
}
