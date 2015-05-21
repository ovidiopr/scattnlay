#ifndef SRC_NMIE_H_
#define SRC_NMIE_H_
//**********************************************************************************//
//    Copyright (C) 2009-2015  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2015  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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

#define VERSION "2.0.0"
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace nmie {
  //Used constants
  const double PI_=3.14159265358979323846;
  // light speed [m s-1]
  const double cc_ = 2.99792458e8;
  // assume non-magnetic (MU=MU0=const) [N A-2]
  const double mu_ = 4.0*PI_*1.0e-7;
  int ScattCoeffs(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const int nmax, std::vector<std::complex<double> >& an, std::vector<std::complex<double> >& bn);
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, const int pl, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nMie(const unsigned int L, std::vector<double>& x, std::vector<std::complex<double> >& m, const unsigned int nTheta, std::vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, std::vector<std::complex<double> >& S1, std::vector<std::complex<double> >& S2);
  int nField(const unsigned int L, const int pl, const std::vector<double>& x, const std::vector<std::complex<double> >& m, const int nmax, const unsigned int ncoord, const std::vector<double>& Xp, const std::vector<double>& Yp, const std::vector<double>& Zp, std::vector<std::vector<std::complex<double> > >& E, std::vector<std::vector<std::complex<double> > >& H);

  class MultiLayerMie {
   public:
    // Run calculation
    void RunMieCalculation();
    void RunFieldCalculation();
    void calcScattCoeffs();

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
    void AddNewLayer(double layer_size, std::complex<double> layer_index);
    // Modify width of the layer
    void SetLayerSize(std::vector<double> layer_size, int layer_position = 0);
    // Modify refractive index of the layer
    void SetLayerIndex(std::vector< std::complex<double> > layer_index, int layer_position = 0);
    // Modify size of all layers
    void SetAllLayersSize(const std::vector<double>& layer_size);
    // Modify refractive index of all layers
    void SetAllLayersIndex(const std::vector< std::complex<double> >& index);
    // Modify scattering (theta) angles
    void SetAngles(const std::vector<double>& angles);
    // Modify coordinates for field calculation
    void SetFieldCoords(const std::vector< std::vector<double> >& coords);
    // Modify PEC layer
    void SetPECLayer(int layer_position = 0);

    // Set a fixed value for the maximun number of terms
    void SetMaxTerms(int nmax);
    // Get maximun number of terms
    int GetMaxTerms() {return nmax_;};

    bool isMieCalculated(){return isMieCalculated_;};
    // Clear layer information
    void ClearLayers();
    void MarkUncalculated();

    // Applied units requests
    double GetSizeParameter();
    double GetLayerWidth(int layer_position = 0);
    std::vector<double> GetLayersSize();
    std::vector<std::complex<double> > GetLayersIndex();
    std::vector<std::array<double, 3> > GetFieldCoords();

    std::vector<std::vector< std::complex<double> > > GetFieldE(){return E_;};   // {X[], Y[], Z[]}
    std::vector<std::vector< std::complex<double> > > GetFieldH(){return H_;};

  protected:
    // Size parameter for all layers
    std::vector<double> size_param_;
    // Refractive index for all layers
    std::vector< std::complex<double> > refractive_index_;
    // Scattering angles for scattering pattern in radians

  private:
    void calcNstop();
    void calcNmax(unsigned int first_layer);

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
    void calcD1D3(std::complex<double> z,
                  std::vector<std::complex<double> >& D1,
                  std::vector<std::complex<double> >& D3);
    void calcPsiZeta(std::complex<double> x,
                     std::vector<std::complex<double> >& Psi,
                     std::vector<std::complex<double> >& Zeta);
    void calcPiTau(const double& costheta,
                   std::vector<double>& Pi, std::vector<double>& Tau);
    void calcSpherHarm(const std::complex<double> Rho, const double Theta, const double Phi,
                       const std::complex<double>& rn, const std::complex<double>& Dn,
                       const double& Pi, const double& Tau, const double& n,
                       std::vector<std::complex<double> >& Mo1n, std::vector<std::complex<double> >& Me1n, 
                       std::vector<std::complex<double> >& No1n, std::vector<std::complex<double> >& Ne1n);
    void calcExpanCoeffs();

    void calcField(const double Rho, const double Theta, const double Phi,
                   std::vector<std::complex<double> >& E, std::vector<std::complex<double> >& H);

    bool isExpCoeffsCalc_ = false;
    bool isScaCoeffsCalc_ = false;
    bool isMieCalculated_ = false;

    std::vector<double> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;

    // with calcNmax(int first_layer);
    int nmax_ = -1;
    int nmax_preset_ = -1;
    // Scattering coefficients
    std::vector<std::complex<double> > an_, bn_;
    std::vector< std::vector<double> > coords_;
    std::vector< std::vector<std::complex<double> > > aln_, bln_, cln_, dln_;
    /// Store result
    double Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    std::vector<std::vector< std::complex<double> > > E_, H_;  // {X[], Y[], Z[]}
    std::vector<std::complex<double> > S1_, S2_;


    //Temporary variables
    std::vector<std::complex<double> > PsiZeta_;


  };  // end of class MultiLayerMie

}  // end of namespace nmie
#endif  // SRC_NMIE_H_
