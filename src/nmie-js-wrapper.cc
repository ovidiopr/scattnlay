//**********************************************************************************//
//    Copyright (C) 2009-2019  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2019  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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
//                                                                                  //
//    @brief  Wrapper to JS
//                                                                                  //
//**********************************************************************************//
#include "nmie-applied.hpp"
#include "nmie-applied-impl.hpp"
#include "nmie-precision.hpp"

using namespace emscripten;

nmie::MultiLayerMieApplied<double> ml_mie;

EMSCRIPTEN_BINDINGS (c) {
        class_<nmie::MultiLayerMieApplied<double>>("nmie")
                .constructor<>()
                .function("SetWavelength", &nmie::MultiLayerMieApplied<double>::SetWavelength)
                .function("AddTargetLayerReIm",&nmie::MultiLayerMieApplied<double>::AddTargetLayerReIm)
                .function("SetModeNmaxAndType",&nmie::MultiLayerMieApplied<double>::SetModeNmaxAndType)
                .function("ClearTarget",&nmie::MultiLayerMieApplied<double>::ClearTarget)
                .function("RunMieCalculation",&nmie::MultiLayerMieApplied<double>::RunMieCalculation)
                .function("RunFieldCalculationPolar",&nmie::MultiLayerMieApplied<double>::RunFieldCalculationPolar)
                .function("GetFieldEabs",&nmie::MultiLayerMieApplied<double>::GetFieldEabs)
                .function("GetQsca",&nmie::MultiLayerMieApplied<double>::GetQsca)
                .function("GetQext",&nmie::MultiLayerMieApplied<double>::GetQext)
                .function("GetQabs",&nmie::MultiLayerMieApplied<double>::GetQabs)
//                .function("bf",&nmie::MultiLayerMieApplied<double>::bf)
        ;
}

//namespace nmie {
  //**********************************************************************************//
  // This function emulates a C call to calculate the actual scattering parameters    //
  // and amplitudes.                                                                  //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nTheta: Number of scattering angles                                            //
  //   Theta: Array containing all the scattering angles where the scattering         //
  //          amplitudes will be calculated                                           //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it              //
  //                                                                                  //
  // Output parameters:                                                               //
  //   Qext: Efficiency factor for extinction                                         //
  //   Qsca: Efficiency factor for scattering                                         //
  //   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
  //   Qbk: Efficiency factor for backscattering                                      //
  //   Qpr: Efficiency factor for the radiation pressure                              //
  //   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
  //   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
  //   S1, S2: Complex scattering amplitudes                                          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
//  int nMieApplied(const unsigned int L, const int pl, std::vector<double> &x, std::vector<std::complex<double> > &m, const unsigned int nTheta, std::vector<double> &Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
//
//    if (x.size() != L || m.size() != L)
//        throw std::invalid_argument("Declared number of layers do not fit x and m!");
//    if (Theta.size() != nTheta)
//        throw std::invalid_argument("Declared number of sample for Theta is not correct!");
//    try {
//      MultiLayerMieApplied<FloatType> ml_mie;
//      ml_mie.SetLayersSize(ConvertVector<FloatType>(x));
//      ml_mie.SetLayersIndex(ConvertComplexVector<FloatType>(m));
//      ml_mie.SetAngles(ConvertVector<FloatType>(Theta));
//      ml_mie.SetPECLayer(pl);
//      ml_mie.SetMaxTerms(nmax);
//
//      ml_mie.RunMieCalculation();
//
//      *Qext = static_cast<double>(ml_mie.GetQext());
//      *Qsca = static_cast<double>(ml_mie.GetQsca());
//      *Qabs = static_cast<double>(ml_mie.GetQabs());
//      *Qbk = static_cast<double>(ml_mie.GetQbk());
//      *Qpr = static_cast<double>(ml_mie.GetQpr());
//      *g = static_cast<double>(ml_mie.GetAsymmetryFactor());
//      *Albedo = static_cast<double>(ml_mie.GetAlbedo());
//      S1 = ConvertComplexVector<double>(ml_mie.GetS1());
//      S2 = ConvertComplexVector<double>(ml_mie.GetS2());
//
//      return ml_mie.GetMaxTerms();
//    } catch(const std::invalid_argument& ia) {
//      // Will catch if  ml_mie fails or other errors.
//      std::cerr << "Invalid argument: " << ia.what() << std::endl;
//      throw std::invalid_argument(ia);
//      return -1;
//    }
//    return 0;
//  }
//  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
//    return nmie::nMieApplied(L, -1, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
//  }
//  int nMieApplied(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
//    return nmie::nMieApplied(L, pl, x, m, nTheta, Theta, -1, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
//  }
//  int nMieApplied(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2) {
//    return nmie::nMieApplied(L, -1, x, m, nTheta, Theta, nmax, Qext, Qsca, Qabs, Qbk, Qpr, g, Albedo, S1, S2);
//  }
  // ********************************************************************** //
  // ********************************************************************** //
  // ********************************************************************** //

//}  // end of namespace nmie
