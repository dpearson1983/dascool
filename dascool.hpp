/**
 * @file dascool.hpp
 * @brief This header file will contain the class structure for David's Cosmology library.
 * 
 * @author David W. Pearson
 * 
 * @date 28 March 2019
 */

#ifndef _DASCOOL_HPP_
#define _DASCOOL_HPP_

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

/**
 * This is used to set the spacing of the frequencies when computing power spectra.
 * @author David W Pearson
 * 
 * @param linear Frequencies are linearly spaced
 * @param logarithmic Frequencies are logarithmically spaced
 * 
 * @date 28 March 2019
 */
enum class spacing{
    linear,
    logarithmic
};

/**
 * This is the class definition.
 * @author David W. Pearson
 * 
 * @param Om_M The density parameter of all matter
 * @param Om_L The density parameter of the cosmological constant
 * @param Om_b The density parameter of baryons
 * @param Om_c The density parameter of cold dark matter
 * @param tau
 * @param T_cmb The cosmic microwave background temperature
 * @param h The Hubble parameter h = H_0/100 km s^{-1} Mpc^{-1}
 * @param n The primordial power spectrum scaling
 * @param sigma_8 The matter fluctuations on scale of 8 Mpc
 * @param k The frequencies at which to calculate power spectra
 * @param acc GSL interpolation acceleration pointer
 * @param r2z GSL spline to correlate comoving distance to redshift
 * 
 * @date 28 March 2019
 */
class dascool{
    double Om_M, Om_L, Om_b, Om_c, tau, T_cmb, h, n, sigma_8; // User provided quantities
    double z_d, z_eq, R_d, R_eq, Theta, k_eq, r_d; // Quantities calculated during initialization
    std::vector<double> k;
    gsl_interp_accel *acc;
    gsl_spline *r2z;
    
    double E(double z);
    
    static double E_inv(double z, void *params);
    
    static double D_1int(double z, void *params);
    
    static double rd_int(double z, void *params);
    
    static double rz(double z, void *params);
    
    void setFrequencies(spacing space, double k_min, double k_max, int num_k);
    
    double calcNoWiggleTransfer(double k_i);
    
    double calcPkPrim(double k_i);
    
    double sigmasqr(double R);
    
    double D_1(double z);
    
    double Growth(double z);
    
    public:
        // Default values from Table 1 of Planck 2018 Results. VI. cosmological parameters
        dascool(double H_O = 67.37, double OmegaM = 0.3147, double OmegaL = 0.6853, double Omegab = 0.049199,
                  double Omegac = 0.26395, double Tau = 0.0540, double TCMB = 2.7255, double ns = 0.9652,
                  double sigma8 = 0.8101);
        
        ~dascool();
        
        double getOmega_M();
        
        double getOmega_L();
        
        double getOmega_bh2();
        
        double getOmega_ch2();
        
        double geth();
        
        double getH_0();
        
        double getz_eq();
        
        double getz_d();
        
        double getr_d();
        
        double R(double z);
        
        double getk_eq();
        
        double H(double z);
        
        double D_A(double z);
        
        double D_V(double z);
        
        double comovingDistance(double z, gsl_integration_workspace *w);
        
        double redshift(double r);
        
        std::vector<std::vector<double>> noWigglePower(spacing space, double k_min, double k_max, int num_k, double z);
        
};
    
#endif
