#pragma once

#include <iostream>
#include <complex>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>		    // Librería para la función gamma logarítmica en los complejos, necesaria para el f_c y el \sigma_l


extern double eta, k;			    // Esta se va a necesitar bastante como parámetro global aquí
std::complex<double> fC_constant;			    // Pre-llamo esta variable para no tener problemas con el compilador

//========================================================================================
//	FUNCIONES GAMMA Y ARGUMENTO DE GAMMA
//========================================================================================
// Como no encontré librerías que me dieran la gamma compleja, tengo que hacer un apaño
// con la gamma logarítmica para poder devoler su valor
std::complex<double> complex_gamma(const std::complex<double>& z) {
    gsl_sf_result lnr, arg;  // Stores ln|Γ(z)| and arg(Γ(z))
    
    // Compute log-Gamma (lnr = log|Γ(z)|, arg = arg(Γ(z)))
    gsl_sf_lngamma_complex_e(z.real(), z.imag(), &lnr, &arg);
    
    // Reconstruct Γ(z) = exp(lnr) * (cos(arg) + i sin(arg))
    double real = std::exp(lnr.val) * std::cos(arg.val);
    double imag = std::exp(lnr.val) * std::sin(arg.val);
    
    return std::complex<double>(real, imag);
}


//  ARGUMENTO DE LA FUNCIÓN GAMMA
//  _____________________________________________________________________
// Función que me devuelve el argumento de la gamma para un complejo z
double argument_complex_gamma(const std::complex<double>& z) {
    gsl_sf_result lnr, arg;  // Stores ln|Γ(z)| and arg(Γ(z))
    
    // Compute log-Gamma (lnr = log|Γ(z)|, arg = arg(Γ(z)))
    gsl_sf_lngamma_complex_e(z.real(), z.imag(), &lnr, &arg);

    return arg.val;
}



//======================================================================
//	CORRIMIENTO Y AMPLITUD DE SCATTERING DE FASE DE COULOMB
//======================================================================

// Corrimiento de fase 
double CoulombPhaseShift(double l){
    return argument_complex_gamma( std::complex<double> (l + 1.0, eta) );
}




// Amplitud de scattering
std::complex<double> CoulombScatteringAmplitude(double theta){
    
    double thetad = theta / 2.0;
    return std::exp(std::complex<double>(0.0, -eta * std::log( std::pow( std::sin(thetad) ,2.0) ))) * 
	   (std::complex<double>( (1.0 / std::pow( std::sin(thetad) ,2.0)) ,0.0)) * 
	   fC_constant;
}





//==================================================================
//	PRUEBA DEL HEADER
//==================================================================
void runHeaderAuxiliar() {
    // Por rendimiento, vamos a obtener la constante que aparece en el corrimiento de fase de Coulomb
    fC_constant = std::complex<double>(-(eta / (2 * k)), 0.0) * std::exp(std::complex<double>(0.0, 1.0) * CoulombPhaseShift(0.0));
    // TODO: Cuando vaya a crear la función Fel, tengo que agregar esta constante en su definición para que se corra y tenerla guardada

    std::complex<double> z (1.0, -1.0);
    std::complex<double> gammaVal = argument_complex_gamma(z);

    std::cout << fC_constant << " || " << CoulombScatteringAmplitude(0.1) << std::endl;
    std::cout << eta << " || " << k << std::endl;
    return;
}
