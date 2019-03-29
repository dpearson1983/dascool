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
 * This function is to be used by GSL to perform the integral needed to compute D_1. This
 * function should not be directly called, it should instead be used to setup a GSL function.
 * @author David W. Pearson
 * 
 * @param z The redshift at which to evaluate the function
 * @param params Pointer to the a struct of intParams
 * 
 * @date 29 March 2019
 */
double dascool::D_1int(double z, void *params) {
    intParams p = *(intParams *)params;
    return (1.0 + z)/sqrt((1.0 + z)*(1.0 + z)*(1.0 + z)*p.OmM + p.OmL);
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
        double pk = T*T*dascool::calcPkPrim(k_i);
        double res = k_i*k_i*k_i*w*w*pk;
        result += w_i[i]*res;
    }
    result *= (k_max - k_min)/2.0;
    return result;
}

/**
 * Calculate the linear growth factor at a given redshift
 * @author David W. Pearson
 * 
 * @param z Redshift at which to calculate
 * 
 * @date 29 March 2019
 */
double dascool::D_1(double z) {
    double intRes, error;
    intParams p;
    p.OmM = dascool::Om_M;
    p.OmL = dascool::Om_L;
    p.Omb = dascool::Om_b;
    p.Tcmb = dascool::T_cmb;
    p.H_0 = 100.0;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000000);
    gsl_function F;
    F.function = &D_1int;
    F.params = &p;
    gsl_integration_qags(&F, z, 1E6, 1E-6, 1E-6, 10000000, w, &intRes, &error);
    gsl_integration_workspace_free(w);
    return dascool::E(z)*intRes;
}

/**
 * Calculate the normalized linear growth factor at a given redshift
 * @author David W. Pearson
 * 
 * @param z Redshift at which to calculate
 * 
 * @date 29 March 2019
 */
double dascool::Growth(double z) {
    return dascool::D_1(z)/dascool::D_1(0);
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
 * @param n The primordial power spectrum slope
 * @param sigma_8 The average fluctuations on scales of 8 Mpc
 * 
 * @date 28 March 2019
 */
dascool::dascool(double H_0, double OmegaM, double OmegaL, double Omegab, double Omegac, double Tau,
                 double TCMB, double ns, double sigma8) {
    dascool::Om_M = OmegaM;
    dascool::Om_L = OmegaL;
    dascool::Om_b = Omegab;
    dascool::Om_c = Omegac;
    dascool::tau = Tau;
    dascool::T_cmb = TCMB;
    dascool::n = ns;
    dascool::sigma_8 = sigma8;
    
    dascool::h = H_0/100.0;
    dascool::Theta = TCMB/2.7;
    
    dascool::acc = gsl_interp_accel_alloc();
    dascool::r2z = gsl_spline_alloc(gsl_interp_cspline, 1001);
    std::vector<double> z;
    std::vector<double> r;
    double dz = 10.0/1000.0;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000000);
    for (int i = 0; i <= 1000; ++i) {
        z.push_back(i*dz);
        r.push_back(dascool::comovingDistance(i*dz, w));
    }
    gsl_spline_init(dascool::r2z, r.data(), z.data(), z.size());
    gsl_integration_workspace_free(w);
}

dascool::~dascool() {
    gsl_spline_free(dascool::r2z);
    gsl_interp_accel_free(dascool::acc);
}

/**
 * Returns the value of Omega_M
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::getOmega_M() {
    return dascool::Om_M;
}

/**
 * Returns the value of Omega_L
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::getOmega_L() {
    return dascool::Om_L;
}

/**
 * Returns the value of Omega_b*h^2
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::getOmega_bh2() {
    return dascool::Om_b*dascool::h*dascool::h;
}

/**
 * Returns the value of Omega_c*h^2
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::getOmega_ch2() {
    return dascool::Om_c*dascool::h*dascool::h;
}

/**
 * Returns the value of h
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::geth() {
    return dascool::h;
}

/**
 * Returns the value of H_0
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */
double dascool::getH_0() {
    return dascool::h*100.0;
}

/**
 * Returns the value of z_eq
 * @author David W. Pearson
 * 
 * @date 29 March 2019
 */
double dascool::getz_eq() {
    return dascool::z_eq;
}

/**
 * Returns the value of z_d
 * @author David W. Pearson
 * 
 * @date 29 March 2019
 */
double dascool::getz_d() {
    return dascool::z_d;
}

/**
 * Returns the value of r_d
 * @author David W. Pearson
 * 
 * @date 29 March 2019
 */
double dascool::getr_d() {
    return dascool::r_d;
}

/**
 * Calculate the ratio of the baryon to photon momentum density at a given redshift
 * @author David W. Pearson
 * 
 * @param z Redshift at which to calculate
 * 
 * @date 29 March 2019
 */
double dascool::R(double z) {
    double Om_bh2 = dascool::Om_b*dascool::h*dascool::h;
    double Theta4 = dascool::Theta*dascool::Theta*dascool::Theta*dascool::Theta;
    return (31.5*Om_bh2)/(Theta4*(z/1000.0));
}

/**
 * Calculate the scale of the particle horizon at z_eq
 * @author David W. Pearson
 * 
 * @param z Redshift at which to calculate.
 * 
 * @date
 */
double dascool::getk_eq() {
//     return sqrt(2.0*dascool::Omega_M()*dascool::H0()*dascool::H0()*dascool::z_eq());
    return (0.0746*dascool::Om_M*dascool::h*dascool::h)/        
                     (dascool::Theta*dascool::Theta);
}

/**
 * Calculate the comoving distance to a given redshift
 * @author David W. Pearson
 * 
 * @param z The redshift to compute the distance to.
 * @param w Point to a GSL integration workspace to remove the overhead of constant allocation and deallocation
 * 
 * @date 29 March 2019
 */
double dascool::comovingDistance(double z, gsl_integration_workspace *w) {
    double D, error;
    intParams p;
    p.OmM = dascool::Om_M;
    p.OmL = dascool::Om_L;
    p.Omb = dascool::Om_b;
    p.Tcmb = dascool::T_cmb;
    p.H_0 = 100.0;
    gsl_function F;
    F.function = &dascool::rz;
    F.params = &p;
    gsl_integration_qags(&F, 0.0, z, 1E-6, 1E-6, 1000000, w, &D, &error);
    return D;
}

/**
 * Calculate the redshift to a give comoving distance
 * @author David W. Pearson
 * 
 * @param r The comoving distance
 * 
 * @date 29 March 2019
 */
double dascool::redshift(double r) {
    return gsl_spline_eval(dascool::r2z, r, dascool::acc);
}

/**
 * Compute linear no-wiggle power spectrum for the defined cosmology with the Eisenstein & Hu 1998 transfer
 * function
 * @author David W. Pearson
 * 
 * @param space The kind of spacing to use for frequencies, k (logarithmic is recommended)
 * @param k_min The minimum frequency
 * @param k_max The maximum frequency
 * @param num_k The total number of frequencies to compute the power spectrum at
 * @param z The redshift to compute the power spectrum at
 * 
 * @date 29 March 2019
 */
std::vector<std::vector<double>> dascool::noWigglePower(spacing space, double k_min, double k_max, int num_k, double z) {
    dascool::setFrequencies(space, k_min, k_max, num_k);
    
    std::vector<std::vector<double>> Pk_NW(2, std::vector<double>(dascool::k.size()));
    double pk_norm = dascool::sigma_8/dascool::sigmasqr(8.0);
    for (int i = 0; i < dascool::k.size(); ++i) {
        double T = calcNoWiggleTransfer(dascool::k[i]);
        double GofZ = dascool::Growth(z);
        Pk_NW[0][i] = dascool::k[i];
        Pk_NW[1][i] = dascool::calcPkPrim(dascool::k[i])*T*T*GofZ*GofZ*pk_norm;
    }
    
    return Pk_NW;
}
