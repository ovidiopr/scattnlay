//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//                                                                                  //
//    This file is part of scattnlay                                                //
//                                                                                  //
//    This program is free software: you can redistribute it and/or modify          //
//    it under the terms of the GNU General Public License as published by          //
//    the Free Software Foundation, either version 3 of the License, or             //
//    (at your option) any later version.                                           //
//                                                                                  //
//    This program is distributed in the hope that it will be useful,               //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                 //
//    GNU General Public License for more details.                                  //
//                                                                                  //
//    The only additional remark is that we expect that all publications            //
//    describing work using this software, or all commercial products               //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#define VERSION "0.3.1"
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace nmie {

  int nMie(int L, std::vector<double>& x, std::vector<std::complex<double> >& m, int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nField(const int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const int ncoord, const std::vector<double>& Xp, const std::vector<double>& Yp, const std::vector<double>& Zp, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H);

  class MultiLayerMie {
   public:
    // Run calculation
    void RunMieCalculations();
    void RunFieldCalculations();

    // Return calculation results
    double GetQext();
    double GetQsca();
    double GetQabs();
    double GetQbk();
    double GetQpr();
    double GetAsymmetryFactor();
    double GetAlbedo();
    std::vector<std::complex<double> > GetS1();
    std::vector<std::complex<double> > GetS2();

    std::vector<std::complex<double> > GetAn(){return an_;};
    std::vector<std::complex<double> > GetBn(){return bn_;}; 

    // Problem definition
    // Add new layer
    void AddNewLayer(double layer_width, std::complex<double> layer_index);
    // Modify width of the layer
    void SetLayerWidth(std::vector<double> layer_width, int layer_position = 0);
    // Modify refractive index of the layer
    void SetLayerIndex(std::vector< std::complex<double> > layer_index, int layer_position = 0);
    // Modify width of all layers
    void SetLayersWidth(std::vector<double> layer_width);
    // Modify refractive index of all layers
    void SetLayerIndex(std::vector< std::complex<double> > layer_index);
    // Set PEC layer
    void SetPECLayer(int layer_position = 0);

    // Set maximun number of terms to be used
    void SetMaxTerms(int nmax);
    // Get maximun number of terms
    int GetMaxTermsUsed() {return nmax_used_;};

    // Clear layer information
    void ClearLayers();

    // Applied units requests
    double GetTotalRadius();
    double GetCoreRadius();
    double GetLayerWidth(int layer_position = 0);
    std::vector<double> GetLayersWidth();
    std::vector<std::complex<double> > GetLayersIndex();  
    std::vector<std::array<double, 3> > GetFieldCoords();

    std::vector<std::vector< std::complex<double> > > GetFieldE(){return E_field_;};   // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<double> > > GetFieldH(){return H_field_;};
  private:
    void Nstop();
    void Nmax(int first_layer);
    void sbesjh(std::complex<double> z, std::vector<std::complex<double> >& jn,
	            std::vector<std::complex<double> >& jnp, std::vector<std::complex<double> >& h1n,
	            std::vector<std::complex<double> >& h1np);
    void sphericalBessel(std::complex<double> z, std::vector<std::complex<double> >& bj,
			             std::vector<std::complex<double> >& by, std::vector<std::complex<double> >& bd);
    std::complex<double> calc_an(int n, double XL, std::complex<double> Ha, std::complex<double> mL,
	                             std::complex<double> PsiXL, std::complex<double> ZetaXL,
				                 std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1);
    std::complex<double> calc_bn(int n, double XL, std::complex<double> Hb, std::complex<double> mL,
	                             std::complex<double> PsiXL, std::complex<double> ZetaXL,
				                 std::complex<double> PsiXLM1, std::complex<double> ZetaXLM1);
    std::complex<double> calc_S1(int n, std::complex<double> an, std::complex<double> bn,
				                 double Pi, double Tau);
    std::complex<double> calc_S2(int n, std::complex<double> an, std::complex<double> bn,
				                 double Pi, double Tau);
    void calcPsiZeta(std::complex<double> x, 
		             std::vector<std::complex<double> > D1,
		             std::vector<std::complex<double> > D3,
		             std::vector<std::complex<double> >& Psi,
		             std::vector<std::complex<double> >& Zeta);
    std::complex<double> calcD1confra(int N, const std::complex<double> z);
    void calcD1D3(std::complex<double> z,
		          std::vector<std::complex<double> >& D1,
		          std::vector<std::complex<double> >& D3);
    void calcSinglePiTau(const double& costheta, std::vector<double>& Pi,
			             std::vector<double>& Tau);
    void calcAllPiTau(std::vector< std::vector<double> >& Pi,
		              std::vector< std::vector<double> >& Tau);
    void ExtScattCoeffs(std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn); 
    void IntScattCoeffs();
    void IntScattCoeffsInit();

    void fieldExt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H);

    void fieldInt(const double Rho, const double Phi, const double Theta, const  std::vector<double>& Pi, const std::vector<double>& Tau, std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H);
    
    bool areIntCoeffsCalc_ = false;
    bool areExtCoeffsCalc_ = false;
    bool isMieCalculated_ = false;
    double wavelength_ = 1.0;
    double total_radius_ = 0.0;

    // Size parameter for all layers
    std::vector<double> width_;
    // Refractive index for all layers
    std::vector< std::complex<double> > index_;
    // Scattering angles for scattering pattern in radians
    std::vector<double> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;

    // with Nmax(int first_layer);
    int nmax_ = -1;
    int nmax_used_ = -1;
    int nmax_preset_ = -1;
    // Scattering coefficients
    std::vector<std::complex<double> > an_, bn_;
    std::vector< std::vector<double> > coords_sp_;
    // TODO: check if l index is reversed will lead to performance
    // boost, if $a^(L+1)_n$ stored in al_n_[n][0], $a^(L)_n$ in
    // al_n_[n][1] and so on...
    // at the moment order is forward!
    std::vector< std::vector<std::complex<double> > > al_n_, bl_n_, cl_n_, dl_n_;
    /// Store result
    double Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    std::vector<std::vector< std::complex<double> > > E_field_, H_field_;  // {X[], Y[], Z[]}
    // Mie efficinecy from each multipole channel.
    std::vector<double> Qsca_ch_, Qext_ch_, Qabs_ch_, Qbk_ch_, Qpr_ch_;
    std::vector<double> Qsca_ch_norm_, Qext_ch_norm_, Qabs_ch_norm_, Qbk_ch_norm_, Qpr_ch_norm_;
    std::vector<std::complex<double> > S1_, S2_;

    //Used constants
    const double PI_=3.14159265358979323846;  
    // light speed [m s-1]
    double const cc_ = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    double const mu_ = 4.0*PI_*1.0e-7;

    //Temporary variables
    std::vector<std::complex<double> > PsiZeta_;


  };  // end of class MultiLayerMie

}  // end of namespace nmie
