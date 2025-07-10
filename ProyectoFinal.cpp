#include <iostream>					// La confiable IOSTREAM
#include <boost/math/special_functions/legendre.hpp>
#include <cmath>					// Para matemáticas en general

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//-------------------------------------------------------
//		CONSTANTES DEL SISTEMA 
//-------------------------------------------------------

#include "headers/CONSTANTS.h"				// Para las constantes del sistema
double E_NN_lab = 85.0;
double M = mN * A;							// Masa de cada núcleo

double gammaL;
void RelativisticTerm(void){
double E_proj = A * E_NN_lab + M;
double gammaL = E_proj / M;
double betaL = std::sqrt(1.0 - (1.0 / (gammaL * gammaL)));
double sqrt_s = std::sqrt(2*M*M + 2*M*(E_proj));
double E_star = sqrt_s / 2.0;
double p_star = std::sqrt(E_star * E_star - M*M);
double K = p_star / hbar;
double v_cm = p_star / E_star;
double v_rel = 2 * v_cm;
double eta = Z*Z*alpha / (v_rel);
double waveNumber = K;
}


// void NonRelativisticTerms(void){
double E_lab = E_NN_lab * A;
double E_cm = E_lab * 0.5;
double K = std::sqrt(2.0 * mu * E_cm ) / hbar;
double v_cm = (K * hbar) / mu;

double eta = (Z * Z * e2) / (hbar * v_cm);
double waveNumber = std::sqrt(2 * mu * E_cm) / hbar;
// }

double k = std::sqrt(2 * mN * E_NN_lab * A) / hbar;
double AB = A * A;
//------------------------------------------------------	
//		CONSTANTES DE LAS FUNCIONES 
//------------------------------------------------------	


#include "headers/AuxiliarFunction.h"			// En donde están las funciones de Coulomb y el bl
double xi_NN_constant = -1.0 / (v_cm * hbar);


// INCLUSIÓN DE LAS LIBRERÍAS PROPIAS
#include "headers/BigIntegral.h"			// En donde se encuentran las integrales feas
#include "headers/BigSummatory.h"			// En este header se define el Sel y el Fel


std::string path = "ExtractedData85MeV.csv";

//------------------------------------------------------	
//		CROSS-SECTION
//------------------------------------------------------	
double CrossSection(double theta, double *params){
	return std::norm( Fel(theta, params) / CoulombScatteringAmplitude(theta)); 
	return std::norm( Fel(theta, params)) / RelRutherfordCrossSection(theta);
}


//------------------------------------------------------	
//		PARA LEER PARÁMETROS DE ENTRADA
//------------------------------------------------------	
std::vector<double> readCSVAngles(void){
	std::ifstream file(path);
	if (!file.is_open()) {
		std::cerr << "Error opening file!" << std::endl;
		return std::vector<double> (1, 0.0);
	}

	std::vector<double> first_column_numbers;
	std::string line;

	while (std::getline(file, line)) {
		std::stringstream ss(line);
		std::string first_field;

		if (std::getline(ss, first_field, ',')) {
		    double value = std::stod(first_field);  // convert string to double
		    first_column_numbers.push_back(value);
		}
	}

	file.close();

	return first_column_numbers;
}



//------------------------------------------------------	
//		FUNCIÓN MAIN
//------------------------------------------------------	
int main(void){
	omp_set_nested(1);
	// double params[4] = {-281.5, 0.0, 0.99, 1};
	double params[4] = {-67.5, 689.6, std::sqrt(0.72), std::sqrt(2.55)}; // E_nn_lab = 85.0


	// double params[5] = {1.0,-67.5, 689.6, 0.72, 2.55}; // E_nn_lab = 85.0
	// for (double j = 1.0; j <= 20; j += 1.0){
	// 	std::cout << "j = " << j << std::endl;
	// 	std::cout << "RealfNN(b =" << j << ") = " << f_NN_real_integrand(j, params) << std::endl;
	// 	std::cout << "ImagfNN(b =" << j << ") = " << f_NN_imaginary_integrand(j, params) << std::endl;
	// }
	
	// return 0;
	
	std::cerr << "eta " << eta << "MeV" << std::endl;
	std::cerr << "K " << K << "MeV" << std::endl;
	std::cerr << "vcm " << v_cm << "MeV" << std::endl;

	std::cerr << "Energía lab. por nucleón: " << E_NN_lab << "MeV" << std::endl;
	// std::cerr << "Energía lab. total: " << E_lab << "MeV" << std::endl;	
	std::cerr << "b Máximo para f_nn: " << bMax << std::endl;
	std::cerr << "q Máximo para Xi: " << qMax << std::endl;
	std::cerr << "Particiones en la integral de Xi: " << N << std::endl;
	std::cerr << "l Máximo a usar: " << lMax << std::endl;



	// for (double i = 0.5; i <= 14; i += 0.25){
	// 	std::cerr << "Calculado theta = " << i << std::endl;
	// 	std::cout << double(i) << "\t" << CrossSection(double(i * pi / 180.0), params) << std::endl;
	// }


		

	std::vector<double> angles = readCSVAngles();

	for (double value : angles){
		std::cerr << "Calculado theta = " << value << std::endl;
		std::cout << value << "\t" << CrossSection(double(value * pi / 180.0), params) << std::endl;
	}


	return 0;
}
