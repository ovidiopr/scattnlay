#ifndef SRC_NMIE_NMIE_H_
#define SRC_NMIE_NMIE_H_
///
/// @file   nmie-wrapper.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Tue Sep  3 00:40:47 2013
/// @copyright 2013 Ladutenko Konstantin
///
/// nmie-wrapper is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// nmie-wrapper is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with nmie-wrapper.  If not, see <http://www.gnu.org/licenses/>.
///
/// nmie-wrapper uses nmie.c from scattnlay by Ovidio Pena
/// <ovidio@bytesfall.com> as a linked library. He has an additional condition to 
/// his library:
//    The only additional condition is that we expect that all publications         //
//    describing  work using this software , or all commercial products             //
//    using it, cite the following reference:                                       //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
///
/// @brief  Wrapper class around nMie function for ease of use
/// 
///
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(EXIT_FAILURE); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif


namespace nmie {

  int nMie_wrapper(int L, const std::vector<double>& x, const std::vector<std::complex<double> >& m, int nTheta, const std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nField(const int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const int ncoord, const std::vector<double>& Xp, const std::vector<double>& Yp, const std::vector<double>& Zp,  std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H)


  class MultiLayerMie {
    // Will throw for any error!
    // SP stands for size parameter units.
   public:
    void GetFailed();
    long iformat = 0;
    bool output = true;
    void prn(double var) {
      do {
	if (!output) break;
	++iformat;
	printf("%23.13e",var);	     
	if (iformat%4 == 0) printf("\n");
      } while (false);
    }
    // Set parameters in applied units 
    void SetWavelength(double wavelength) {wavelength_ = wavelength;};
    // It is possible to set only a multilayer target to run calculaitons.
    // For many runs it can be convenient to separate target and coating layers.
    // Per layer
    void AddTargetLayer(double layer_width, std::complex<double> layer_index);
    void AddCoatingLayer(double layer_width, std::complex<double> layer_index);
    // For all layers
    void SetTargetWidth(std::vector<double> width);
    void SetTargetIndex(std::vector< std::complex<double> > index);
    void SetTargetPEC(double radius);
    void SetCoatingWidth(std::vector<double> width);
    void SetCoatingIndex(std::vector< std::complex<double> > index);
    void SetFieldPoints(std::vector< std::array<double,3> > coords);

    //Set parameters in size parameter units
    void SetWidthSP(const std::vector<double>& width);
    void SetIndexSP(const std::vector< std::complex<double> >& index);
    void SetFieldPointsSP(const std::vector< std::vector<double> >& coords_sp);

    // Set common parameters
    void SetAnglesForPattern(double from_angle, double to_angle, int samples);
    void SetAngles(const std::vector<double>& angles);
    std::vector<double> GetAngles();
    void SetPEC(int layer_position = 0);  // By default set PEC layer to be the first one
    
    void SetMaxTermsNumber(int nmax);
    int GetMaxTermsUsed() {return nmax_used_;};
    
    void ClearTarget();
    void ClearCoating();
    void ClearLayers();
    void ClearAllDesign(); //Layers + SP + index_

    // Applied units requests
    double GetTotalRadius();
    double GetTargetRadius();
    double GetCoatingWidth();
    std::vector<double>                  GetTargetLayersWidth();
    std::vector< std::complex<double> >  GetTargetLayersIndex();
    std::vector<double>                  GetCoatingLayersWidth();
    std::vector< std::complex<double> >  GetCoatingLayersIndex();
    std::vector< std::array<double,3> >   GetFieldPoints();
    std::vector<std::vector< std::complex<double> > >  GetFieldE();
    std::vector<std::vector< std::complex<double> > >  GetFieldH();
    std::vector< std::vector<double> >   GetSpectra(double from_WL, double to_WL,
                                                   int samples);  // ext, sca, abs, bk
    double GetRCSext();
    double GetRCSsca();
    double GetRCSabs();
    double GetRCSbk();
    std::vector<double> GetPatternEk();
    std::vector<double> GetPatternHk();
    std::vector<double> GetPatternUnpolarized();
    


    // Size parameter units
    std::vector<double>                  GetLayerWidthSP();
    // Same as to get target and coating index
    std::vector< std::complex<double> >  GetLayerIndex();  
    std::vector< std::array<double,3> >   GetFieldPointsSP();
    // Do we need normalize field to size parameter?
    /* std::vector<std::vector<std::complex<double> > >  GetFieldESP(); */
    /* std::vector<std::vector<std::complex<double> > >  GetFieldHSP(); */
    std::vector< std::array<double,5> >   GetSpectraSP(double from_SP, double to_SP,
						       int samples);  // WL,ext, sca, abs, bk
    double GetQext();
    double GetQsca();
    double GetQabs();
    double GetQbk();
    double GetQpr();
    std::vector<double> GetQsca_channel();
    std::vector<double> GetQabs_channel();
    std::vector<double> GetQsca_channel_normalized();
    std::vector<double> GetQabs_channel_normalized();

    double GetAsymmetryFactor();
    double GetAlbedo();
    std::vector<std::complex<double> > GetS1();
    std::vector<std::complex<double> > GetS2();
    std::vector<double> GetPatternEkSP();
    std::vector<double> GetPatternHkSP();
    std::vector<double> GetPatternUnpolarizedSP();
    
    // Run calculation
    void RunMieCalculations();
    void RunFieldCalculations();

    // Output results (data file + python script to plot it with matplotlib)
    void PlotSpectra();
    void PlotSpectraSP();
    void PlotField();
    void PlotFieldSP();
    void PlotPattern();
    void PlotPatternSP();

  private:
    void ConvertToSP();
    void GenerateSizeParameter();
    void GenerateIndex();
    void InitMieCalculations();

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
    void calcPsiZeta(double x, 
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
    void calcAllPiTau( std::vector< std::vector<double> >& Pi,
		    std::vector< std::vector<double> >& Tau);
    void ScattCoeffs(std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn); 
    void fieldExt( double Rho, double Phi, double Theta, std::vector<double> Pi, std::vector<double> Tau,
		  std::vector<std::complex<double> > an, std::vector<std::complex<double> > bn,
		  std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H);
    
    bool isMieCalculated_ = false;
    double wavelength_ = 1.0;
    double total_radius_ = 0.0;
    /// Width and index for each layer of the structure
    std::vector<double> target_width_, coating_width_;
    std::vector< std::complex<double> > target_index_, coating_index_;
    /// Size parameters for all layers
    std::vector<double> size_parameter_;
    /// Complex index values for each layers.
    std::vector< std::complex<double> > index_;
    /// Scattering angles for RCS pattern in radians
    std::vector<double> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;
    // Set nmax_ manualy with SetMaxTermsNumber(int nmax) or in ScattCoeffs(..)
    // with Nmax(int first_layer);
    int nmax_ = -1;
    int nmax_used_ = -1;
    int nmax_preset_ = -1;
    // Scattering coefficients
    std::vector<std::complex<double> > an_, bn_;
    std::vector< std::vector<double> > coords_sp_;
    /// Store result
    double Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    // Mie efficinecy from each multipole channel.
    std::vector<double> Qsca_ch_, Qext_ch_, Qabs_ch_, Qbk_ch_, Qpr_ch_;
    std::vector<double> Qsca_ch_norm_, Qext_ch_norm_, Qabs_ch_norm_, Qbk_ch_norm_, Qpr_ch_norm_;
    std::vector<std::complex<double> > S1_, S2_;

    //Used constants
    const double PI=3.14159265358979323846;  
    // light speed [m s-1]
    double const cc = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    double const mu = 4.0*PI*1.0e-7;

    //Temporary variables
    std::vector<std::complex<double> > PsiZeta_;


  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_NMIE_H_
