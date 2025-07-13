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
double E_NN_lab = 200.0;
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
std::vector<double> readCSVAngles(std::string path){
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
int main(int argc, char *argv[]){

	omp_set_nested(1);
	// double params[4] = {-281.5, 0.0, 0.99, 1};
	//
	
	bool fit = 1;
	std::string path;
	double maxAngle;

	double params[4];
	if (argc == 1){
		params[0] = -67.5; params[1] =  689.6; params[2] = std::sqrt(0.72); params[3] = std::sqrt(2.55); // E_nn_lab = 85.0
		path = "ExtractedData30MeV.csv";
	}  
	else if (argc == 5){
		params[0] = std::stod(argv[1]); params[1] = std::stod(argv[2]); params[2] = std::stod(argv[3]); params[3] = std::stod(argv[4]);
		path = "ExtractedData200MeV.csv";
	}
	else if (argc == 7){
		params[0] = std::stod(argv[3]); params[1] = std::stod(argv[4]); params[2] = std::stod(argv[5]); params[3] = std::stod(argv[6]);
		if (std::string(argv[1]) == "-f" || std::string(argv[2]) == "--fit"){
			// std::string path = "ExtractedData85MeV.csv";
			path = argv[2];
		}
		else if(std::string(argv[1]) == "-g" || std::string(argv[2]) == "--graph"){
			fit = 0;
			maxAngle = std::stod(argv[2]);
		}
	}
	else {
		std::cerr << "Brother, use the proper params wtf" << std::endl;
		return 1;
	}

	
	if (fit == 1){
		std::vector<double> angles = readCSVAngles(path);

		for (double value : angles){
			// std::cerr << "Calculado theta = " << value << std::endl;
			std::cout << value << "\t" << CrossSection(double(value * pi / 180.0), params) << std::endl;
		}
	} else {
		for (double value = 0.5; value <= maxAngle; value += 0.05){
			std::cout << value << "\t" << CrossSection(double(value * pi / 180.0), params) << std::endl;
		}
	}


	return 0;
	
}
