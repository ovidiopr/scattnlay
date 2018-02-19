#ifndef SRC_NMIE_APPLIED_HPP_
#define SRC_NMIE_APPLIED_HPP_
//**********************************************************************************//
//    Copyright (C) 2009-2018  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2018  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//    using it, cite at least one of the following references:                      //
//    [1] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [2] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
//    You should have received a copy of the GNU General Public License             //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.         //
//**********************************************************************************//

#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "nmie.hpp"
#include "nmie-impl.hpp"


namespace nmie {

  int nMieApplied(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMieApplied(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);


  template <typename FloatType = double>
  class MultiLayerMieApplied : public MultiLayerMie<FloatType> {
    // Will throw for any error!
   public:
    void RunMieCalculation();
    void GetFailed();
    long iformat = 0;
    bool output = true;
    void prn(FloatType var) {
      do {
	if (!output) break;
	++iformat;
	printf("%23.13e",var);	     
	if (iformat%4 == 0) printf("\n");
      } while (false);
    }
    // Set parameters in applied units 
    void SetWavelength(FloatType wavelength) {wavelength_ = wavelength;};
    // It is possible to set only a multilayer target to run calculaitons.
    // For many runs it can be convenient to separate target and coating layers.
    // Per layer
    void AddTargetLayer(FloatType layer_width, std::complex<FloatType> layer_index);
    void AddCoatingLayer(FloatType layer_width, std::complex<FloatType> layer_index);
    // For all layers
    void SetTargetWidth(std::vector<FloatType> width);
    void SetTargetIndex(std::vector< std::complex<FloatType> > index);
    void SetTargetPEC(FloatType radius);
    void SetCoatingWidth(std::vector<FloatType> width);
    void SetCoatingIndex(std::vector< std::complex<FloatType> > index);
    void SetFieldPoints(std::vector< std::array<FloatType,3> > coords);

    //Set parameters in size parameter units
    void SetWidthSP(const std::vector<FloatType>& width);
    void SetIndexSP(const std::vector< std::complex<FloatType> >& index);
    void SetFieldPointsSP(const std::vector< std::vector<FloatType> >& coords_sp);

    // Set common parameters
    void SetAnglesForPattern(FloatType from_angle, FloatType to_angle, int samples);
    std::vector<FloatType> GetAngles();
    
    void ClearTarget();
    void ClearCoating();
    void ClearLayers();
    void ClearAllDesign(); //Layers + SP + index_

    // Applied units requests
    FloatType GetTotalRadius();
    FloatType GetTargetRadius();
    FloatType GetCoatingWidth();
    std::vector<FloatType>                  GetTargetLayersWidth();
    std::vector< std::complex<FloatType> >  GetTargetLayersIndex();
    std::vector<FloatType>                  GetCoatingLayersWidth();
    std::vector< std::complex<FloatType> >  GetCoatingLayersIndex();
    std::vector< std::vector<FloatType> >   GetFieldPoints();
    std::vector< std::vector<FloatType> > GetSpectra(FloatType from_WL, FloatType to_WL, int samples);  // ext, sca, abs, bk
    FloatType GetRCSext();
    FloatType GetRCSsca();
    FloatType GetRCSabs();
    FloatType GetRCSbk();
    std::vector<FloatType> GetPatternEk();
    std::vector<FloatType> GetPatternHk();
    std::vector<FloatType> GetPatternUnpolarized();

    // Size parameter units
    std::vector<FloatType> GetLayerWidthSP();
    // Same as to get target and coating index
    std::vector< std::complex<FloatType> > GetLayerIndex();  
    std::vector< std::array<FloatType,3> > GetFieldPointsSP();
    // Do we need normalize field to size parameter?
    /* std::vector<std::vector<std::complex<FloatType> > >  GetFieldESP(); */
    /* std::vector<std::vector<std::complex<FloatType> > >  GetFieldHSP(); */
    std::vector< std::array<FloatType,5> > GetSpectraSP(FloatType from_SP, FloatType to_SP, int samples);  // WL,ext, sca, abs, bk


    std::vector<FloatType> GetPatternEkSP();
    std::vector<FloatType> GetPatternHkSP();
    std::vector<FloatType> GetPatternUnpolarizedSP();

    void GetExpanCoeffs
      (std::vector< std::vector<std::complex<FloatType> > >& aln,
       std::vector< std::vector<std::complex<FloatType> > >& bln,
       std::vector< std::vector<std::complex<FloatType> > >& cln,
       std::vector< std::vector<std::complex<FloatType> > >& dln);


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

    void sbesjh(std::complex<FloatType> z, std::vector<std::complex<FloatType> >& jn,
	            std::vector<std::complex<FloatType> >& jnp, std::vector<std::complex<FloatType> >& h1n,
	            std::vector<std::complex<FloatType> >& h1np);
    void sphericalBessel(std::complex<FloatType> z, std::vector<std::complex<FloatType> >& bj,
			             std::vector<std::complex<FloatType> >& by, std::vector<std::complex<FloatType> >& bd);
    std::complex<FloatType> calcD1confra(int N, const std::complex<FloatType> z);
    
    FloatType wavelength_ = 1.0;
    FloatType total_radius_ = 0.0;
    /// Width and index for each layer of the structure
    std::vector<FloatType> target_width_, coating_width_;
    std::vector< std::complex<FloatType> > target_index_, coating_index_;

    std::vector< std::vector<FloatType> > coords_sp_;




  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_APPLIED_HPP
