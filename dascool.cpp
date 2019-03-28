/**
 * @file dascool.cpp
 * @brief This implementation file will contain the class function definitions for David's Cosmology library.
 * 
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */

#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "dascool.hpp"

#define c 299792.458 // Speed of light in km s^{-1}
#define G 6.6740831E-11 // Newton's gravitational constant
#define sigma 5.67037321E-8
#define convert 1.0502650413E-30
#define pi 3.14159265359 // It's pi

const double w_i[] = {0.096540088514728, 0.096540088514728, 0.095638720079275, 0.095638720079275,
                     0.093844399080805, 0.093844399080805, 0.091173878695764, 0.091173878695764,
                     0.087652093004404, 0.087652093004404, 0.083311924226947, 0.083311924226947,
                     0.078193895787070, 0.078193895787070, 0.072345794108849, 0.072345794108849,
                     0.065822222776362, 0.065822222776362, 0.058684093478536, 0.058684093478536,
                     0.050998059262376, 0.050998059262376, 0.042835898022227, 0.042835898022227,
                     0.034273862913021, 0.034273862913021, 0.025392065309262, 0.025392065309262,
                     0.016274394730906, 0.016274394730906, 0.007018610009470, 0.007018610009470};

const double x_i[] = {-0.048307665687738, 0.048307665687738, -0.144471961582796, 0.144471961582796,
                     -0.239287362252137, 0.239287362252137, -0.331868602282128, 0.331868602282128,
                     -0.421351276130635, 0.421351276130635, -0.506899908932229, 0.506899908932229,
                     -0.587715757240762, 0.587715757240762, -0.663044266930215, 0.663044266930215,
                     -0.732182118740290, 0.732182118740290, -0.794483795967942, 0.794483795967942,
                     -0.849367613732570, 0.849367613732570, -0.896321155766052, 0.896321155766052,
                     -0.934906075937739, 0.934906075937739, -0.964762255587506, 0.964762255587506,
                     -0.985611511545268, 0.985611511545268, -0.997263861849481, 0.997263861849481};

/**
 * This struct is used to pass parameters to functions that GSL will integrate
 * @author David W. Pearson
 * 
 * @param OmM The density parameter of matter
 * @param OmL The density parameter of the cosmological constant
 * @param Omb The density parameter of baryons
 * @param Tcmb The cosmic microwave background temperature
 * @param H_0 The Hubble Constant
 * 
 * @date 28 March 2019
 */
struct intParams{
    double OmM, OmL, Omb, Tcmb, H_0;
};

/**
 * This function returns value of the normalized Hubble function at redshift z in a flat LCDM universe
 * @author David W. Pearson
 * 
 * @param z The redshift at which to evaluate the function
 * 
 * @date 28 March 2019
 */
double dascool::E(double z) {
    return sqrt((1.0 + z)*(1.0 + z)*(1.0 + z)*dascool::Om_M + dascool::Om_L);
}

/**
 * This function is to be used by GSL to integrate the inverse of the normalized Hubble function. This
 * function should not be directly called, it should instead be used to setup a GSL function.
 * @author David W. Pearson
 * 
 * @param z The redshift at which to evaluate the function
 * @param params Pointer to the a struct of intParams
 * 
 * @date 28 March 2019
 */
double dascool::E_inv(double z, void *params) {
    intParams p = *(intParams *)params;
    return 1.0/sqrt((1.0 + z)*(1.0 + z)*(1.0 + z)*p.OmM + p.OmL);
}

/**
 * This function is to be used by GSL to integrate in order to calculate the drag radius. This function
 * should not be directly called, it should instead be used to setup a GSL function.
 * @author David W. Pearson
 * 
 * @param z The redshift at which to evaluate the function
 * @param params Pointer to the a struct of intParams
 * 
 * @date 28 March 2019
 */
double dascool::rd_int(double z, void *params) {
    intParams p = *(intParams *)params;
    double coeff2 = (9.0*c*c*c*p.H_0*p.H_0*p.Omb*convert)/(128.0*sigma*pi*G*p.Tcmb*p.Tcmb*p.Tcmb*p.Tcmb);
    double E = sqrt((1.0 + z)*(1.0 + z)*(1.0 + z)*p.OmM + p.OmL);
    double term1 = (coeff2/((1.0 + z))) + 1.0;
    return 1.0/(sqrt(term1)*E);
}

/**
 * This function is to be used by GSL to integrate in order to calculate the comoving distance. This function
 * should not be directly called, it should instead be used to setup a GSL function.
 * @author David W. Pearson
 * 
 * @param z The redshift at which to evaluate the function
 * @param params Pointer to the a struct of intParams
 * 
 * @date 28 March 2019
 */
double dascool::rz(double z, void *params) {
    intParams p = *(intParams *)params;
    double D = c/(100.0*sqrt(p.OmM*(1.0 + z)*(1.0 + z)*(1.0 + z) + p.OmL));
    return D;
}

/**
 * This function will initialize the frequency vector class member.
 * @author David W. Pearson
 * 
 * @param space The kind of spacing to use for the frequencies
 * @param k_min The minimum frequency
 * @param k_max The maximum frequency
 * @param num_k The number of frequencies
 * 
 * @date 28 March 2019
 */
void dascool::setFrequencies(spacing space, double k_min, double k_max, int num_k) {
    if (space == spacing::linear) {
        double dk = (k_max - k_min)/(num_k - 1.0);
        for (int i = 0; i < num_k; ++i) {
            double k_i = k_min + i*dk;
            dascool::k.push_back(k_i);
        }
    } else if (space == spacing::logarithmic) {
        double dk = (log10(k_max) - log10(k_min))/(num_k - 1.0);
        for (int i = 0; i < num_k; ++i) {
            double k_i = log10(k_min) + i*dk;
            k_i = std::pow(10, k_i);
            dascool::k.push_back(k_i);
        }
    }
}

/**
 * This function will initialize the no-wiggle (i.e. no BAO features) Transfer function from Eisenstein & Hu 
 * 1998
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::calcNoWiggleTransfer(double k_i) {
    double w_m = dascool::Om_M*dascool::h*dascool::h;
    double w_b = dascool::Om_b*dascool::h*dascool::h;
    double fb = dascool::Om_b/dascool::Om_M;
    double fc = (dascool::Om_M - dascool::Om_b)/dascool::Om_M;
    double z_eq = 2.5E4*w_m/(dascool::Theta*dascool::Theta*dascool::Theta*dascool::Theta);
    double k_eq = 7.46E-2*w_m/(dascool::Theta*dascool::Theta);
    double b1 = 0.313*pow(w_m, -0.419)*(1.0 + 0.607*pow(w_m, 0.674));
    double b2 = 0.238*pow(w_m, 0.223);
    double z_d = 1291.0*pow(w_m, 0.251)*(1.0 + b1*pow(w_b, b2))/(1.0 + 0.659*pow(w_m, 0.828));
    double R_d = (31.5*w_b/(dascool::Theta*dascool::Theta*dascool::Theta*dascool::Theta))*(1000.0/z_d);
    double R_eq = (31.5*w_b/(dascool::Theta*dascool::Theta*dascool::Theta*dascool::Theta))*(1000.0/z_eq);
    double sh_d = (2.0/(3.0*k_eq))*sqrt(6.0/R_eq)*log((sqrt(1.0 + R_d) + sqrt(R_d + R_eq))/(1.0 + sqrt(R_eq)));
    double alpha_gamma = 1.0 - 0.328*log(431.0*w_m)*fb + 0.38*log(22.3*w_m)*fb*fb;
    
    double quadTerm = (0.43*k_i*sh_d)*(0.43*k_i*sh_d)*(0.43*k_i*sh_d)*(0.43*k_i*sh_d);
    double gamma_eff = dascool::Om_M*dascool::h*(alpha_gamma + (1.0 - alpha_gamma)/(1.0 + quadTerm));
    double q = k_i*dascool::Theta*dascool::Theta/gamma_eff;
    double L = log(2.0*exp(1.0) + 1.8*q);
    double C = 14.2 + 731.0/(1.0 + 62.5*q);
    return L/(L + C*q*q);
}

/**
 * This function will calculate the primordial power spectrum
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::calcPkPrim(double k_i) {
        return std::pow(k_i, dascool::n);
}

/**
 * This function will calculate the matter fluctuations within a sphere of radius R
 * @author David W. Pearson
 * 
 * @param R The radius in which to calculate fluctuations
 * 
 * @date 28 March 2019
 */
double dascool::sigmasqr(double R) {
    // Integrate using 32 Gaussian quadrature, which is probably overkill
    double k_min = dascool::k[0];
    double k_max = dascool::k[dascool::k.size() - 1];
    double result = 0.0;
    for (int i = 0; i < 32; ++i) {
        double k_i = ((k_max - k_min)/2.0)*x_i[i] + (k_max + k_min)/2.0;
        double x = k_i*R;
        double w = 3.0*(std::sin(x) - x*std::cos(x))/(x*x*x);
        double T = dascool::calcNoWiggleTransfer(k_i);
        double pk = T*T*dascool::pk_prim[i];
        double res = k_i*k_i*k_i*w*w*pk;
        result += w_i[i]*res
    }
    result *= (k_max - k_min)/2.0;
    return result;
}

/**
 * This is the class initializer. Defaults to Planck 2018 LCDM cosmology.
 * @author David W. Pearson
 * 
 * @param H_0 The Hubble constant
 * @param OmegaM The density parameter of all matter
 * @param OmegaL The density parameter of the cosmological constant
 * @param Omegab The density parameter of baryons
 * @param Omegac The density parameter of cold dark matter
 * @param Tau The reionization optical depth
 * @param TCMB The cosmic microwave background temperature
 * 
