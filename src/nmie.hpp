#ifndef SRC_NMIE_HPP_
#define SRC_NMIE_HPP_
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

#define VERSION "2.2"
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <boost/math/constants/constants.hpp>
namespace nmie {
  int ScattCoeffs(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nmax, std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn);
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nField(const unsigned int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const unsigned int ncoord, const std::vector<double>& Xp, const std::vector<double>& Yp, const std::vector<double>& Zp, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H);

  template <typename FloatType = double>
  class MultiLayerMie {    
   public:
    //Used constants TODO! Change to boost PI
    const double PI_=3.14159265358979323846;
    // light speed [m s-1]
    const double cc_ = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    const double mu_ = 4.0*PI_*1.0e-7;
    // Run calculation
    void RunMieCalculation();
    void RunFieldCalculation();
    void calcScattCoeffs();

    // Return calculation results
    FloatType GetQext();
    FloatType GetQsca();
    FloatType GetQabs();
    FloatType GetQbk();
    FloatType GetQpr();
    FloatType GetAsymmetryFactor();
    FloatType GetAlbedo();
    std::vector<std::complex<FloatType> > GetS1();
    std::vector<std::complex<FloatType> > GetS2();

    std::vector<std::complex<FloatType> > GetAn(){return an_;};
    std::vector<std::complex<FloatType> > GetBn(){return bn_;};

    // Problem definition
    // Modify size of all layers
    void SetLayersSize(const std::vector<FloatType>& layer_size);
    // Modify refractive index of all layers
    void SetLayersIndex(const std::vector< std::complex<FloatType> >& index);
    // Modify scattering (theta) angles
    void SetAngles(const std::vector<FloatType>& angles);
    // Modify coordinates for field calculation
    void SetFieldCoords(const std::vector< std::vector<FloatType> >& coords);
    // Modify index of PEC layer
    void SetPECLayer(int layer_position = 0);

    // Set a fixed value for the maximun number of terms
    void SetMaxTerms(int nmax);
    // Get maximun number of terms
    int GetMaxTerms() {return nmax_;};

    bool isMieCalculated(){return isMieCalculated_;};
    // Clear layer information
    void ClearLayers();
    void MarkUncalculated();

    // Read parameters
    // Get total size parameter of particle
    FloatType GetSizeParameter();
    // Returns size of all layers
    std::vector<FloatType> GetLayersSize(){return size_param_;};
    // Returns refractive index of all layers
    std::vector<std::complex<FloatType> > GetLayersIndex(){return refractive_index_;};
    // Returns scattering (theta) angles
    std::vector<FloatType> GetAngles(){return theta_;};
    // Returns coordinates used for field calculation
    std::vector<std::vector<FloatType> > GetFieldCoords(){return coords_;};
    // Returns index of PEC layer
    int GetPECLayer(){return PEC_layer_position_;};

    std::vector<std::vector< std::complex<FloatType> > > GetFieldE(){return E_;};   // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<FloatType> > > GetFieldH(){return H_;};
    // Get fields in spherical coordinates.
    std::vector<std::vector< std::complex<FloatType> > > GetFieldEs(){return E_;};   // {rho[], teha[], phi[]}
    std::vector<std::vector< std::complex<FloatType> > > GetFieldHs(){return H_;};

  protected:
    // Size parameter for all layers
    std::vector<FloatType> size_param_;
    // Refractive index for all layers
    std::vector< std::complex<FloatType> > refractive_index_;
    // Scattering coefficients
    std::vector<std::complex<FloatType> > an_, bn_;
    std::vector< std::vector<std::complex<FloatType> > > aln_, bln_, cln_, dln_;
    void calcExpanCoeffs();
    // Points for field evaluation
    std::vector< std::vector<FloatType> > coords_;

  private:
    void calcNstop();
    void calcNmax(unsigned int first_layer);

    std::complex<FloatType> calc_an(int n, FloatType XL, std::complex<FloatType> Ha, std::complex<FloatType> mL,
                                 std::complex<FloatType> PsiXL, std::complex<FloatType> ZetaXL,
                                 std::complex<FloatType> PsiXLM1, std::complex<FloatType> ZetaXLM1);
    std::complex<FloatType> calc_bn(int n, FloatType XL, std::complex<FloatType> Hb, std::complex<FloatType> mL,
                                 std::complex<FloatType> PsiXL, std::complex<FloatType> ZetaXL,
                                 std::complex<FloatType> PsiXLM1, std::complex<FloatType> ZetaXLM1);
    std::complex<FloatType> calc_S1(int n, std::complex<FloatType> an, std::complex<FloatType> bn,
                                 FloatType Pi, FloatType Tau);
    std::complex<FloatType> calc_S2(int n, std::complex<FloatType> an, std::complex<FloatType> bn,
                                 FloatType Pi, FloatType Tau);
    void calcD1D3(std::complex<FloatType> z,
                  std::vector<std::complex<FloatType> >& D1,
                  std::vector<std::complex<FloatType> >& D3);
    void calcPsiZeta(std::complex<FloatType> x,
                     std::vector<std::complex<FloatType> >& Psi,
                     std::vector<std::complex<FloatType> >& Zeta);
    void calcPiTau(const FloatType& costheta,
                   std::vector<FloatType>& Pi, std::vector<FloatType>& Tau);
    void calcSpherHarm(const std::complex<FloatType> Rho, const FloatType Theta, const FloatType Phi,
                       const std::complex<FloatType>& rn, const std::complex<FloatType>& Dn,
                       const FloatType& Pi, const FloatType& Tau, const FloatType& n,
                       std::vector<std::complex<FloatType> >& Mo1n, std::vector<std::complex<FloatType> >& Me1n, 
                       std::vector<std::complex<FloatType> >& No1n, std::vector<std::complex<FloatType> >& Ne1n);

    void calcField(const FloatType Rho, const FloatType Theta, const FloatType Phi,
                   std::vector<std::complex<FloatType> >& E, std::vector<std::complex<FloatType> >& H);

    bool isExpCoeffsCalc_ = false;
    bool isScaCoeffsCalc_ = false;
    bool isMieCalculated_ = false;

    // Scattering angles for scattering pattern in radians
    std::vector<FloatType> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;

    // with calcNmax(int first_layer);
    int nmax_ = -1;
    int nmax_preset_ = -1;
    /// Store result
    FloatType Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    std::vector<std::vector< std::complex<FloatType> > > E_, H_;  // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<FloatType> > > Es_, Hs_;  // {X[], Y[], Z[]}
    std::vector<std::complex<FloatType> > S1_, S2_;


    //Temporary variables
    std::vector<std::complex<FloatType> > PsiZeta_;


  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_HPP_
