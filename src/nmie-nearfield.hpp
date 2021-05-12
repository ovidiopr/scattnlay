#ifndef SRC_NMIE_NEARFIELD_HPP_
#define SRC_NMIE_NEARFIELD_HPP_
//**********************************************************************************//
//    Copyright (C) 2009-2021  Ovidio Pena <ovidio@bytesfall.com>                   //
//    Copyright (C) 2013-2021  Konstantin Ladutenko <kostyfisik@gmail.com>          //
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

//**********************************************************************************//
// This class implements the algorithm for a multilayered sphere described by:      //
//    [1] W. Yang, "Improved recursive algorithm for light scattering by a          //
//        multilayered sphere,” Applied Optics, vol. 42, Mar. 2003, pp. 1710-1720.  //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//    [2] O. Pena and U. Pal, "Scattering of electromagnetic radiation by           //
//        a multilayered sphere," Computer Physics Communications,                  //
//        vol. 180, Nov. 2009, pp. 2348-2354.                                       //
//    [3] K. Ladutenko, U. Pal, A. Rivera, and O. Pena-Rodriguez, "Mie              //
//        calculation of electromagnetic near-field for a multilayered              //
//        sphere," Computer Physics Communications, vol. 214, May 2017,             //
//        pp. 225-230.                                                              //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>

//#include "nmie.hpp"

namespace nmie {
  //class implementation

  //**********************************************************************************//
  // This function calculates the expansion coefficients inside the particle,         //
  // required to calculate the near-field parameters.                                 //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send -1                          //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to -1 and the function will calculate it.             //
  //                                                                                  //
  // Output parameters:                                                               //
  //   aln, bln, cln, dln: Complex scattering amplitudes inside the particle          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  template <typename FloatType>
  void MultiLayerMie<FloatType>::calcExpanCoeffs() {
    if (!isScaCoeffsCalc_)
      throw std::invalid_argument("(calcExpanCoeffs) You should calculate external coefficients first!");

    isExpCoeffsCalc_ = false;

    std::complex<FloatType> c_one(1.0, 0.0), c_zero(0.0, 0.0);

    const int L = refractive_index_.size();

    aln_.resize(L + 1);
    bln_.resize(L + 1);
    cln_.resize(L + 1);
    dln_.resize(L + 1);
    for (int l = 0; l <= L; l++) {
      aln_[l].resize(nmax_);
      bln_[l].resize(nmax_);
      cln_[l].resize(nmax_);
      dln_[l].resize(nmax_);
    }

    // Yang, paragraph under eq. A3
    // a^(L + 1)_n = a_n, d^(L + 1) = 1 ...
    for (int n = 0; n < nmax_; n++) {
      aln_[L][n] = an_[n];
      bln_[L][n] = bn_[n];
      cln_[L][n] = c_one;
      dln_[L][n] = c_one;
    }

    std::vector<std::complex<FloatType> > D1z(nmax_ + 1), D1z1(nmax_ + 1), D3z(nmax_ + 1), D3z1(nmax_ + 1);
    std::vector<std::complex<FloatType> > Psiz(nmax_ + 1), Psiz1(nmax_ + 1), Zetaz(nmax_ + 1), Zetaz1(nmax_ + 1);
    std::complex<FloatType> denomZeta, denomPsi, T1, T2, T3, T4;

    auto &m = refractive_index_;
    std::vector< std::complex<FloatType> > m1(L);

    for (int l = 0; l < L - 1; l++) m1[l] = m[l + 1];
    m1[L - 1] = std::complex<FloatType> (1.0, 0.0);

    std::complex<FloatType> z, z1;
    for (int l = L - 1; l >= 0; l--) {
      if (l <= PEC_layer_position_) { // We are inside a PEC. All coefficients must be zero!!!
        for (int n = 0; n < nmax_; n++) {
          // aln
          aln_[l][n] = c_zero;
          // bln
          bln_[l][n] = c_zero;
          // cln
          cln_[l][n] = c_zero;
          // dln
          dln_[l][n] = c_zero;
        }
      } else { // Regular material, just do the calculation
        z = size_param_[l]*m[l];
        z1 = size_param_[l]*m1[l];

        calcD1D3(z, D1z, D3z);
        calcD1D3(z1, D1z1, D3z1);
        calcPsiZeta(z, Psiz, Zetaz);
        calcPsiZeta(z1, Psiz1, Zetaz1);

        for (int n = 0; n < nmax_; n++) {
          int n1 = n + 1;

          denomZeta = Zetaz[n1]*(D1z[n1] - D3z[n1]);
          denomPsi  =  Psiz[n1]*(D1z[n1] - D3z[n1]);

          T1 =  aln_[l + 1][n]*Zetaz1[n1] - dln_[l + 1][n]*Psiz1[n1];
          T2 = (bln_[l + 1][n]*Zetaz1[n1] - cln_[l + 1][n]*Psiz1[n1])*m[l]/m1[l];

          T3 = (dln_[l + 1][n]*D1z1[n1]*Psiz1[n1] - aln_[l + 1][n]*D3z1[n1]*Zetaz1[n1])*m[l]/m1[l];
          T4 =  cln_[l + 1][n]*D1z1[n1]*Psiz1[n1] - bln_[l + 1][n]*D3z1[n1]*Zetaz1[n1];

          // aln
          aln_[l][n] = (D1z[n1]*T1 + T3)/denomZeta;
          // bln
          bln_[l][n] = (D1z[n1]*T2 + T4)/denomZeta;
          // cln
          cln_[l][n] = (D3z[n1]*T2 + T4)/denomPsi;
          // dln
          dln_[l][n] = (D3z[n1]*T1 + T3)/denomPsi;
        }  // end of all n
      }  // end PEC condition
    }  // end of all l

    // Check the result and change  aln_[0][n] and aln_[0][n] for exact zero
    for (int n = 0; n < nmax_; ++n) {
      if (cabs(aln_[0][n]) < 1e-10) aln_[0][n] = 0.0;
      else {
        //throw std::invalid_argument("Unstable calculation of aln_[0][n]!");
        std::cout<< std::setprecision(100)
                 << "Warning: Potentially unstable calculation of aln[0]["
                 << n << "] = "<< aln_[0][n] <<std::endl;
        aln_[0][n] = 0.0;
      }
      if (cabs(bln_[0][n]) < 1e-10) bln_[0][n] = 0.0;
      else {
        //throw std::invalid_argument("Unstable calculation of bln_[0][n]!");
        std::cout<< std::setprecision(100)
                 << "Warning: Potentially unstable calculation of bln[0]["
                 << n << "] = "<< bln_[0][n] <<std::endl;
        bln_[0][n] = 0.0;
      }
    }

    isExpCoeffsCalc_ = true;
  }  // end of   void MultiLayerMie::calcExpanCoeffs()


  template <typename FloatType>
  void MultiLayerMie<FloatType>::convertFieldsFromSphericalToCartesian() {
    long total_points = coords_polar_.size();
    E_.clear(); H_.clear();
    Eabs_.clear(); Habs_.clear();
    for (int point=0; point < total_points; point++) {
      auto Theta = coords_polar_[point][1];
      auto Phi = coords_polar_[point][2];
      auto Es = Es_[point];
      auto Hs = Hs_[point];
      using nmm::sin;
      using nmm::cos;
      E_.push_back({ sin(Theta)*cos(Phi)*Es[0] + cos(Theta)*cos(Phi)*Es[1] - sin(Phi)*Es[2],
                     sin(Theta)*sin(Phi)*Es[0] + cos(Theta)*sin(Phi)*Es[1] + cos(Phi)*Es[2],
                     cos(Theta)*Es[0] - sin(Theta)*Es[1]});
      H_.push_back({ sin(Theta)*cos(Phi)*Hs[0] + cos(Theta)*cos(Phi)*Hs[1] - sin(Phi)*Hs[2],
                     sin(Theta)*sin(Phi)*Hs[0] + cos(Theta)*sin(Phi)*Hs[1] + cos(Phi)*Hs[2],
                     cos(Theta)*Hs[0] - sin(Theta)*Hs[1]});
      Eabs_.push_back(vabs(E_.back()));
      Habs_.push_back(vabs(H_.back()));
    }

  }  // end of void MultiLayerMie::convertFieldsFromSphericalToCartesian()
  //**********************************************************************************//
  // This function calculates the electric (E) and magnetic (H) fields inside and     //
  // around the particle.                                                             //
  //
  // Main troubles of near-field evaluations originate from special functions
  // evaluation, so we expect that nmax needed for the convergence is the size
  // of Psi vector.
  //                                                                                  //
  // Input parameters (coordinates of the point):                                     //
  //   Rho: Radial distance                                                           //
  //   Phi: Azimuthal angle                                                           //
  //   Theta: Polar angle                                                             //
  //   mode_n: mode order.                                                            //
  //          -1 - use all modes (all_)                                               //
  //           1 - use dipole mode only                                               //
  //           2 - use quadrupole mode only                                           //
  //           ...                                                                    //
  //   mode_type: only used when mode_n != -1                                         //
  //          0 - electric only                                                       //
  //          1 - magnetic only                                                       //
  //                                                                                  //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic fields                                     //
  //**********************************************************************************//
  template <typename FloatType>  template <typename evalType>
  void MultiLayerMie<FloatType>::calcFieldByComponents(const evalType Rho,
                                const evalType Theta, const evalType Phi,
                                const std::vector<std::complex<evalType> > &Psi,
                                const std::vector<std::complex<evalType> > &D1n,
                                const std::vector<std::complex<evalType> > &Zeta,
                                const std::vector<std::complex<evalType> > &D3n,
                                const std::vector<evalType> &Pi,
                                const std::vector<evalType> &Tau,
                                std::vector<std::complex<evalType> > &E,
                                std::vector<std::complex<evalType> > &H)  {
    auto nmax = Psi.size() - 1;
    std::complex<evalType> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
    // Vector containing precomputed integer powers of i to avoid computation
    std::vector<std::complex<evalType> > ipow = {c_one, c_i, -c_one, -c_i};
    std::vector<std::complex<evalType> > M3o1n(3), M3e1n(3), N3o1n(3), N3e1n(3);
    std::vector<std::complex<evalType> > M1o1n(3), M1e1n(3), N1o1n(3), N1e1n(3);

    std::complex<evalType> ml;

    // Initialize E and H
    for (int i = 0; i < 3; i++) {
      E[i] = c_zero;
      H[i] = c_zero;
    }

    unsigned int l;
    GetIndexAtRadius(Rho, ml, l);

    for (unsigned int n = 0; n < nmax; n++) {
      int n1 = n + 1;
      auto rn = static_cast<evalType>(n1);

      // using BH 4.12 and 4.50
      calcSpherHarm(Rho*ml, Theta, Phi, Psi[n1], D1n[n1], Pi[n], Tau[n], rn, M1o1n, M1e1n, N1o1n, N1e1n);
      calcSpherHarm(Rho*ml, Theta, Phi, Zeta[n1], D3n[n1], Pi[n], Tau[n], rn, M3o1n, M3e1n, N3o1n, N3e1n);

      // Total field in the lth layer: eqs. (1) and (2) in Yang, Appl. Opt., 42 (2003) 1710-1720
      std::complex<evalType> En = ipow[n1 % 4]
      *static_cast<evalType>((rn + rn + 1.0)/(rn*rn + rn));
      std::complex<evalType> Ediff, Hdiff;
      std::complex<FloatType> Ediff_ft, Hdiff_ft;
      for (int i = 0; i < 3; i++) {
        std::complex<evalType> aln = ConvertComplex<evalType>(aln_[l][n]);
        std::complex<evalType> bln = ConvertComplex<evalType>(bln_[l][n]);
        std::complex<evalType> cln = ConvertComplex<evalType>(cln_[l][n]);
        std::complex<evalType> dln = ConvertComplex<evalType>(dln_[l][n]);
        Ediff = En*(      cln*M1o1n[i] - c_i*dln*N1e1n[i]
                         + c_i*aln*N3e1n[i] -     bln*M3o1n[i]);
        Hdiff = En*(     -dln*M1e1n[i] - c_i*cln*N1o1n[i]
                         + c_i*bln*N3o1n[i] +     aln*M3e1n[i]);
        Ediff_ft = ConvertComplex<FloatType>(Ediff);
        Hdiff_ft = ConvertComplex<FloatType>(Hdiff);
        if ( nmm::isnan(Ediff_ft.real()) || nmm::isnan(Ediff_ft.imag()) ||
             nmm::isnan(Hdiff_ft.real()) || nmm::isnan(Hdiff_ft.imag()) ) {
          std::cout << "Unexpected truncation during near-field evaluation at n = "<< n
                    << " (of total nmax = "<<nmax<<")!!!"<<std::endl;
          break;
        }
        if (mode_n_ == Modes::kAll) {
          // electric field E [V m - 1] = EF*E0
          E[i] += Ediff;
          H[i] += Hdiff;
          continue;
        }
        if (n1 == mode_n_) {
          if (mode_type_ == Modes::kElectric || mode_type_ == Modes::kAll) {
            E[i] += En*( -c_i*dln*N1e1n[i]
                        + c_i*aln*N3e1n[i]);

            H[i] += En*(-dln*M1e1n[i]
                        +aln*M3e1n[i]);
            //std::cout << mode_n_;
          }
          if (mode_type_ == Modes::kMagnetic  || mode_type_ == Modes::kAll) {
            E[i] += En*(  cln*M1o1n[i]
                        - bln*M3o1n[i]);

            H[i] += En*( -c_i*cln*N1o1n[i]
                        + c_i*bln*N3o1n[i]);
            //std::cout << mode_n_;
          }
          //std::cout << std::endl;
        }
        //throw std::invalid_argument("Error! Unexpected mode for field evaluation!\n mode_n="+std::to_string(mode_n)+", mode_type="+std::to_string(mode_type)+"\n=====*****=====");
      }
      if (nmm::isnan(Ediff_ft.real()) || nmm::isnan(Ediff_ft.imag()) ||
          nmm::isnan(Hdiff_ft.real()) || nmm::isnan(Hdiff_ft.imag())
          ) break;
    }  // end of for all n

    // magnetic field
    std::complex<evalType> hffact = ml/static_cast<evalType>(cc_*mu_);
    for (int i = 0; i < 3; i++) {
      H[i] = hffact*H[i];
    }
   }  // end of MultiLayerMie::calcFieldByComponents(...)


  //**********************************************************************************//
  // This function calculates complex electric and magnetic field in the surroundings //
  // and inside the particle.                                                         //
  //                                                                                  //
  // Input parameters:                                                                //
  //   L: Number of layers                                                            //
  //   pl: Index of PEC layer. If there is none just send 0 (zero)                    //
  //   x: Array containing the size parameters of the layers [0..L-1]                 //
  //   m: Array containing the relative refractive indexes of the layers [0..L-1]     //
  //   nmax: Maximum number of multipolar expansion terms to be used for the          //
  //         calculations. Only use it if you know what you are doing, otherwise      //
  //         set this parameter to 0 (zero) and the function will calculate it.       //
  //   ncoord: Number of coordinate points                                            //
  //   Coords: Array containing all coordinates where the complex electric and        //
  //           magnetic fields will be calculated                                     //
  //   mode_n: mode order.                                                            //
  //          -1 - use all modes (all_)                                               //
  //           1 - use dipole mode only                                               //
  //           2 - use quadrupole mode only                                           //
  //           ...                                                                    //
  //   mode_type: only used when mode_n != -1                                         //
  //          0 - electric only                                                       //
  //          1 - magnetic only                                                       //
  //                                                                                  //
  // Output parameters:                                                               //
  //   E, H: Complex electric and magnetic field at the provided coordinates          //
  //                                                                                  //
  // Return value:                                                                    //
  //   Number of multipolar expansion terms used for the calculations                 //
  //**********************************************************************************//
  template <typename FloatType>
  void MultiLayerMie<FloatType>::RunFieldCalculation() {
    // Calculate scattering coefficients an_ and bn_
    calcScattCoeffs();
    // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
    calcExpanCoeffs();

    Es_.clear(); Hs_.clear(); coords_polar_.clear();
    long total_points = coords_[0].size();
    for (int point = 0; point < total_points; point++) {
      const FloatType &Xp = coords_[0][point];
      const FloatType &Yp = coords_[1][point];
      const FloatType &Zp = coords_[2][point];

      // Convert to spherical coordinates
      auto Rho = nmm::sqrt(pow2(Xp) + pow2(Yp) + pow2(Zp));
      // If Rho=0 then Theta is undefined. Just set it to zero to avoid problems
      auto Theta = (Rho > 0.0) ? nmm::acos(Zp/Rho) : 0.0;
      // std::atan2 should take care of any special cases, e.g.  Xp=Yp=0, etc.
      auto Phi = nmm::atan2(Yp,Xp);
      coords_polar_.push_back({Rho, Theta, Phi});
      // Avoid convergence problems due to Rho too small
      if (Rho < 1e-5) Rho = 1e-5;

      //*******************************************************//
      // external scattering field = incident + scattered      //
      // BH p.92 (4.37), 94 (4.45), 95 (4.50)                  //
      // assume: medium is non-absorbing; refim = 0; Uabs = 0  //
      //*******************************************************//

      // This array contains the fields in spherical coordinates
      std::vector<std::complex<FloatType> > Es(3), Hs(3);

      // Do the actual calculation of electric and magnetic field
      std::vector<std::complex<FloatType> > Psi(nmax_ + 1), D1n(nmax_ + 1), Zeta(nmax_ + 1), D3n(nmax_ + 1);
      std::vector<FloatType> Pi(nmax_), Tau(nmax_);
      std::complex<FloatType> ml;
      GetIndexAtRadius(Rho, ml);

      // Calculate logarithmic derivative of the Ricatti-Bessel functions
      calcD1D3(Rho*ml, D1n, D3n);
      // Calculate Ricatti-Bessel functions
      calcPsiZeta(Rho*ml, Psi, Zeta);
      // Calculate angular functions Pi and Tau
      calcPiTau(nmm::cos(Theta), Pi, Tau);

      calcFieldByComponents(Rho, Theta, Phi, Psi, D1n, Zeta, D3n, Pi, Tau, Es, Hs);
      Es_.push_back(Es);
      Hs_.push_back(Hs);
    }  // end of for all field coordinates
    convertFieldsFromSphericalToCartesian();
  }  //  end of MultiLayerMie::RunFieldCalculation()

template <typename FloatType>
double eval_delta(const unsigned int steps, const double from_value, const double to_value) {
  auto delta = std::abs(from_value - to_value);
  if (steps < 2) return delta;
  delta /= static_cast<double>(steps-1);
  // We have a limited double precision evaluation of special functions, typically it is 1e-10.
  if ( (2.*delta)/std::abs(from_value+to_value) < 1e-9)
    throw std::invalid_argument("Error! The step is too fine, not supported!");
  return delta;
}


// ml - refractive index
// l - Layer number
template <typename FloatType> template <typename evalType>
void MultiLayerMie<FloatType>::GetIndexAtRadius(const evalType Rho,
                                                std::complex<evalType> &ml,
                                                unsigned int &l) {
  l = 0;
  if (Rho > size_param_.back()) {
    l = size_param_.size();
    ml = std::complex<evalType>(1.0, 0.0);
  } else {
    for (int i = size_param_.size() - 1; i >= 0 ; i--) {
      if (Rho <= size_param_[i]) {
        l = i;
      }
    }
    ml = ConvertComplex<evalType>(refractive_index_[l]);
  }
}
template <typename FloatType> template <typename evalType>
void MultiLayerMie<FloatType>::GetIndexAtRadius(const evalType Rho,
                                                std::complex<evalType> &ml) {
  unsigned int l;
  GetIndexAtRadius(Rho, ml, l);
}

template <typename FloatType>
void MultiLayerMie<FloatType>::calcMieSeriesNeededToConverge(const FloatType Rho) {
  auto required_near_field_nmax = calcNmax(Rho);
  SetMaxTerms(required_near_field_nmax);
  // Calculate scattering coefficients an_ and bn_
  calcScattCoeffs();
  // We might be limited with available machine precision
  available_maximal_nmax_ = nmax_;
  // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
  calcExpanCoeffs();
}


template <typename FloatType>
void MultiLayerMie<FloatType>::calcRadialOnlyDependantFunctions(const double from_Rho, const double to_Rho,
                                                                const bool isIgnoreAvailableNmax,
                                                                std::vector<std::vector<std::complex<FloatType> > > &Psi,
                                                                std::vector<std::vector<std::complex<FloatType> > > &D1n,
                                                                std::vector<std::vector<std::complex<FloatType> > > &Zeta,
                                                                std::vector<std::vector<std::complex<FloatType> > > &D3n) {
  auto radius_points = Psi.size();
  std::vector<std::vector<std::complex<FloatType> > > PsiZeta(radius_points);
  double delta_Rho = eval_delta<double>(radius_points, from_Rho, to_Rho);
  for (unsigned int j=0; j < radius_points; j++) {
    auto Rho = static_cast<FloatType>(from_Rho + j*delta_Rho);
//    if (Rho < 1e-5) Rho = 1e-5; // TODO do we need this?.
    int near_field_nmax = calcNmax(Rho);
    // Skip if not enough terms in Mie series (i.e. required near field nmax > available terms )
    if (near_field_nmax > available_maximal_nmax_ && !isIgnoreAvailableNmax) continue;
    if (near_field_nmax > available_maximal_nmax_)  near_field_nmax = available_maximal_nmax_;
    Psi[j].resize(near_field_nmax + 1); D1n[j].resize(near_field_nmax + 1);
    Zeta[j].resize(near_field_nmax + 1); D3n[j].resize(near_field_nmax + 1);
    PsiZeta[j].resize(near_field_nmax + 1);
    std::complex<FloatType> ml;
    GetIndexAtRadius(Rho, ml);
    auto z = Rho*ml;
    evalDownwardD1(z, D1n[j]);
    evalUpwardPsi(z,  D1n[j], Psi[j]);
    evalUpwardD3 (z, D1n[j], D3n[j], PsiZeta[j]);
    for (unsigned int k = 0; k < Zeta[j].size(); k++) {
      Zeta[j][k] = PsiZeta[j][k]/Psi[j][k];
    }
  }

}


// input parameters:
//         outer_arc_points: will be increased to the nearest power of 2.
template <typename FloatType>
void MultiLayerMie<FloatType>::RunFieldCalculationPolar(const int outer_arc_points,
                                                        const int radius_points,
                                                        const double from_Rho, const double to_Rho,
                                                        const double from_Theta, const double to_Theta,
                                                        const double from_Phi, const double to_Phi,
                                                        const bool isIgnoreAvailableNmax) {
  if (from_Rho > to_Rho || from_Theta > to_Theta || from_Phi > to_Phi
      || outer_arc_points < 1 || radius_points < 1
      || from_Rho < 0.)
    throw std::invalid_argument("Error! Invalid argument for RunFieldCalculationPolar() !");
  auto nmax_old = nmax_;
  int theta_points = 0, phi_points = 0;
  if (to_Theta-from_Theta > to_Phi-from_Phi) {
    theta_points = outer_arc_points;
    phi_points =  static_cast<int>((to_Phi-from_Phi)/(to_Theta-from_Theta) * outer_arc_points);
  } else {
    phi_points = outer_arc_points;
    theta_points =  static_cast<int>((to_Theta-from_Theta)/(to_Phi-from_Phi) * outer_arc_points);
  }
  if (theta_points == 0) theta_points = 1;
  if (phi_points == 0) phi_points = 1;
  calcMieSeriesNeededToConverge(to_Rho);

  std::vector<std::vector<FloatType> >  Pi(theta_points), Tau(theta_points);
  calcPiTauAllTheta(from_Theta, to_Theta, Pi, Tau);

  std::vector<std::vector<std::complex<FloatType> > > Psi(radius_points), D1n(radius_points),
      Zeta(radius_points), D3n(radius_points), PsiZeta(radius_points);
  calcRadialOnlyDependantFunctions(from_Rho, to_Rho, isIgnoreAvailableNmax,
                                   Psi, D1n, Zeta, D3n);

  double delta_Rho = eval_delta<double>(radius_points, from_Rho, to_Rho);
  double delta_Theta = eval_delta<double>(theta_points, from_Theta, to_Theta);
  double delta_Phi = eval_delta<double>(phi_points, from_Phi, to_Phi);
  Es_.clear(); Hs_.clear(); coords_polar_.clear();
  for (int j=0; j < radius_points; j++) {
    auto Rho = from_Rho + j * delta_Rho;
    std::vector< std::complex<double> > Psi_dp = ConvertComplexVector<double>(Psi[j]);
    std::vector< std::complex<double> > Zeta_dp = ConvertComplexVector<double>(Zeta[j]);
    std::vector< std::complex<double> > D1n_dp = ConvertComplexVector<double>(D1n[j]);
    std::vector< std::complex<double> > D3n_dp = ConvertComplexVector<double>(D3n[j]);
    for (int i = 0; i < theta_points; i++) {
      auto Theta = from_Theta + i * delta_Theta;
      std::vector<double> Pi_dp = ConvertVector<double>(Pi[i]);
      std::vector<double> Tau_dp = ConvertVector<double>(Tau[i]);
      for (int k = 0; k < phi_points; k++) {

        auto Phi = from_Phi + k * delta_Phi;
        coords_polar_.push_back({Rho, Theta, Phi});

        // This array contains the fields in spherical coordinates
//        std::vector<std::complex<FloatType> > Es(3), Hs(3);
        std::vector<std::complex<double> > Es(3), Hs(3);
        calcFieldByComponents(
            Rho, Theta, Phi,
                              Psi_dp, D1n_dp, Zeta_dp, D3n_dp,
                              Pi_dp, Tau_dp, Es, Hs
//        Rho, Theta, Phi,
//            Psi[j], D1n[j], Zeta[j], D3n[j],
//            Pi[i], Tau[i], Es, Hs
        );
        Es_.push_back(ConvertComplexVector<FloatType>(Es));
        Hs_.push_back(ConvertComplexVector<FloatType>(Hs));
      }
    }
  }
  convertFieldsFromSphericalToCartesian();
  nmax_ = nmax_old;
}

// Python interface


}  // end of namespace nmie
#endif  // SRC_NMIE_NEARFIELD_HPP_
