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

#define VERSION "0.3.1"
#include <array>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

namespace nmie {
  int ScattCoeffs(const unsigned int L, const int pl, vector<double>& x, vector<complex<double> >& m, const int nmax, vector<complex<double> >& an, vector<complex<double> >& bn);
  int nMie(const unsigned int L, const int pl, vector<double>& x, vector<complex<double> >& m, const unsigned int nTheta, vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2);
  int nMie(const unsigned int L, vector<double>& x, vector<complex<double> >& m, const unsigned int nTheta, vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2);
  int nMie(const unsigned int L, const int pl, vector<double>& x, vector<complex<double> >& m, const unsigned int nTheta, vector<double>& Theta, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2);
  int nMie(const unsigned int L, vector<double>& x, vector<complex<double> >& m, const unsigned int nTheta, vector<double>& Theta, const int nmax, double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, vector<complex<double> >& S1, vector<complex<double> >& S2);
  int nField(const unsigned int L, const int pl, const vector<double>& x, const vector<complex<double> >& m, const int nmax, const unsigned int ncoord, const vector<double>& Xp, const vector<double>& Yp, const vector<double>& Zp, vector<vector<complex<double> > >& E, vector<vector<complex<double> > >& H);

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
    vector<complex<double> > GetS1();
    vector<complex<double> > GetS2();

    vector<complex<double> > GetAn(){return an_;};
    vector<complex<double> > GetBn(){return bn_;};

    // Problem definition
    // Add new layer
    void AddNewLayer(double layer_size, complex<double> layer_index);
    // Modify width of the layer
    void SetLayerSize(vector<double> layer_size, int layer_position = 0);
    // Modify refractive index of the layer
    void SetLayerIndex(vector< complex<double> > layer_index, int layer_position = 0);
    // Modify size of all layers
    void SetAllLayersSize(const vector<double>& layer_size);
    // Modify refractive index of all layers
    void SetAllLayersIndex(const vector< complex<double> >& index);
    // Modify scattering (theta) angles
    void SetAngles(const vector<double>& angles);
    // Modify coordinates for field calculation
    void SetFieldCoords(const vector< vector<double> >& coords);
    // Modify PEC layer
    void SetPECLayer(int layer_position = 0);

    // Set a fixed value for the maximun number of terms
    void SetMaxTerms(int nmax);
    // Get maximun number of terms
    int GetMaxTerms() {return nmax_;};

    // Clear layer information
    void ClearLayers();

    // Applied units requests
    double GetSizeParameter();
    double GetLayerWidth(int layer_position = 0);
    vector<double> GetLayersSize();
    vector<complex<double> > GetLayersIndex();
    vector<array<double, 3> > GetFieldCoords();

    vector<vector< complex<double> > > GetFieldE(){return E_;};   // {X[], Y[], Z[]}
    vector<vector< complex<double> > > GetFieldH(){return H_;};
  private:
    void calcNstop();
    void calcNmax(unsigned int first_layer);

    complex<double> calc_an(int n, double XL, complex<double> Ha, complex<double> mL,
                            complex<double> PsiXL, complex<double> ZetaXL,
                            complex<double> PsiXLM1, complex<double> ZetaXLM1);
    complex<double> calc_bn(int n, double XL, complex<double> Hb, complex<double> mL,
                            complex<double> PsiXL, complex<double> ZetaXL,
                            complex<double> PsiXLM1, complex<double> ZetaXLM1);
    complex<double> calc_S1(int n, complex<double> an, complex<double> bn,
                            double Pi, double Tau);
    complex<double> calc_S2(int n, complex<double> an, complex<double> bn,
                            double Pi, double Tau);
    void calcD1D3(complex<double> z,
                  vector<complex<double> >& D1,
                  vector<complex<double> >& D3);
    void calcPsiZeta(complex<double> x,
                     vector<complex<double> >& Psi,
                     vector<complex<double> >& Zeta);
    void calcPiTau(const double& costheta,
                   vector<double>& Pi, vector<double>& Tau);
    void calcSpherHarm(const complex<double> Rho, const double Theta, const double Phi,
                       const complex<double>& rn, const complex<double>& Dn,
                       const double& Pi, const double& Tau, const double& n,
                       vector<complex<double> >& Mo1n, vector<complex<double> >& Me1n, 
                       vector<complex<double> >& No1n, vector<complex<double> >& Ne1n);
    void calcExpanCoeffs();

    void calcField(const double Rho, const double Theta, const double Phi,
                   vector<complex<double> >& E, vector<complex<double> >& H);

    bool isExpCoeffsCalc_ = false;
    bool isScaCoeffsCalc_ = false;
    bool isMieCalculated_ = false;

    // Size parameter for all layers
    vector<double> size_param_;
    // Refractive index for all layers
    vector< complex<double> > refractive_index_;
    // Scattering angles for scattering pattern in radians
    vector<double> theta_;
    // Should be -1 if there is no PEC.
    int PEC_layer_position_ = -1;

    // with calcNmax(int first_layer);
    int nmax_ = -1;
    int nmax_preset_ = -1;
    // Scattering coefficients
    vector<complex<double> > an_, bn_;
    vector< vector<double> > coords_;
    // TODO: check if l index is reversed will lead to performance
    // boost, if $a^(L+1)_n$ stored in aln_[n][0], $a^(L)_n$ in
    // aln_[n][1] and so on...
    // at the moment order is forward!
    vector< vector<complex<double> > > aln_, bln_, cln_, dln_;
    /// Store result
    double Qsca_ = 0.0, Qext_ = 0.0, Qabs_ = 0.0, Qbk_ = 0.0, Qpr_ = 0.0, asymmetry_factor_ = 0.0, albedo_ = 0.0;
    vector<vector< complex<double> > > E_, H_;  // {X[], Y[], Z[]}
    vector<complex<double> > S1_, S2_;

    //Used constants
    const double PI_=3.14159265358979323846;
    // light speed [m s-1]
    double const cc_ = 2.99792458e8;
    // assume non-magnetic (MU=MU0=const) [N A-2]
    double const mu_ = 4.0*PI_*1.0e-7;

    //Temporary variables
    vector<complex<double> > PsiZeta_;


  };  // end of class MultiLayerMie

}  // end of namespace nmie
