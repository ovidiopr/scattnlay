#ifNdef SRC_NMIE_NMIE_WRAPPER_H_
#define SRC_NMIE_NMIE_WRAPPER_H_
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
  class MultiLayerMie {
    // Will throw for any error!
    // SP stands for size parameter units.
   public:
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
    void SetCoatingWidth(std::vector<double> width);
    void SetCoatingIndex(std::vector< std::complex<double> > index);
    void SetFieldPoints(std::vector< std::array<double,3> > coords);

    //Set parameters in size parameter units
    void SetWidthSP(std::vector<double> width);
    void SetIndexSP(std::vector< std::complex<double> > index);
    void SetFieldPointsSP(std::vector< std::array<double,3> > coords);

    // Set common parameters
    void SetAnglesForPattern(double from_angle, double to_angle, int samples);
    std::vector<double> GetAngles();
    
    void ClearTarget();
    void ClearCoating();
    void ClearLayers();

    // Applied units requests
    double GetTotalRadius();
    double GetTargetRadius();
    double GetCoatingWidth();
    std::vector<double>                  GetTargetLayersWidth();
    std::vector< std::complex<double> >  GetTargetLayersIndex();
    std::vector<double>                  GetCoatingLayersWidth();
    std::vector< std::complex<double> >  GetCoatingLayersIndex();
    std::vector< std::vector<double> >   GetFieldPoints();
    std::vector<std::array< std::complex<double>,3 > >  GetFieldE();
    std::vector<std::array< std::complex<double>,3 > >  GetFieldH();
    std::vector< std::array<double,4> >   GetSpectra(double from_WL, double to_WL,
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
    std::vector< std::vector<double> >   GetFieldPointsSP();
    // Do we need normalize field to size parameter?
    /* std::vector<std::vector<std::complex<double> > >  GetFieldESP(); */
    /* std::vector<std::vector<std::complex<double> > >  GetFieldHSP(); */
    std::vector< std::array<double,4> >   GetSpectraSP(double from_SP, double to_SP,
						       int samples);  // ext, sca, abs, bk
    double GetQext();
    double GetQsca();
    double GetQabs();
    double GetQbk();
    double GetQpr();
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

  private:
    const double PI=3.14159265358979323846;
    void GenerateSizeParameter();
    void GenerateIndex();
    double wavelength_ = 1.0;
    double total_radius_ = 0.0;
    /// Width and index for each layer of the structure
    std::vector<double> target_width_, coating_width_;
    std::vector< std::complex<double> > target_index_, coating_index_;
    /// Size parameters for all layers
    std::vector<double> size_parameter_;
    /// Complex index values for each layers.
    std::vector< std::complex<double> > index_;
  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_NMIE_WRAPPER_H_
