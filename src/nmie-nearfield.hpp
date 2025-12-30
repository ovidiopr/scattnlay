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
#include "special-functions-impl.hpp"

namespace nmie {
template <typename FloatType, typename Engine>
struct NearFieldKernelHelper {
  using RealV = typename Engine::RealV;
  using ComplexV = typename Engine::ComplexV;

  static inline void calcSpherHarm(
      const ComplexV Rho,
      const RealV Theta,
      const RealV Phi,
      const ComplexV rn,
      const ComplexV Dn,
      const RealV Pi,
      const RealV Tau,
      const RealV n,
      ComplexV Mo1n[3],
      ComplexV Me1n[3],
      ComplexV No1n[3],
      ComplexV Ne1n[3]) {
    
    auto zero = Engine::set(0.0);
    auto c_zero = Engine::make_complex(zero, zero);

    auto sin_Phi = Engine::sin(Phi);
    auto cos_Phi = Engine::cos(Phi);
    auto sin_Theta = Engine::sin(Theta);

    // Mo1n
    Mo1n[0] = c_zero;
    // Mo1n[1] = cos_Phi * Pi * rn / Rho;
    auto term1 = Engine::div(Engine::mul(Engine::make_complex(Engine::mul(cos_Phi, Pi), zero), rn), Rho);
    Mo1n[1] = term1;
    // Mo1n[2] = -sin_Phi * Tau * rn / Rho;
    auto term2 = Engine::div(Engine::mul(Engine::make_complex(Engine::mul(Engine::sub(zero, sin_Phi), Tau), zero), rn), Rho);
    Mo1n[2] = term2;

    // Me1n
    Me1n[0] = c_zero;
    // Me1n[1] = -sin_Phi * Pi * rn / Rho;
    Me1n[1] = Engine::div(Engine::mul(Engine::make_complex(Engine::mul(Engine::sub(zero, sin_Phi), Pi), zero), rn), Rho);
    // Me1n[2] = -cos_Phi * Tau * rn / Rho;
    Me1n[2] = Engine::div(Engine::mul(Engine::make_complex(Engine::mul(Engine::sub(zero, cos_Phi), Tau), zero), rn), Rho);

    // No1n
    // No1n[0] = sin_Phi * (n * n + n) * sin_Theta * Pi * rn / Rho / Rho;
    auto n_sq_plus_n = Engine::add(Engine::mul(n, n), n);
    auto factor = Engine::mul(sin_Phi, Engine::mul(n_sq_plus_n, Engine::mul(sin_Theta, Pi)));
    No1n[0] = Engine::div(Engine::div(Engine::mul(Engine::make_complex(factor, zero), rn), Rho), Rho);
    
    // No1n[1] = sin_Phi * Tau * Dn * rn / Rho;
    No1n[1] = Engine::div(Engine::mul(Engine::mul(Engine::make_complex(Engine::mul(sin_Phi, Tau), zero), Dn), rn), Rho);
    
    // No1n[2] = cos_Phi * Pi * Dn * rn / Rho;
    No1n[2] = Engine::div(Engine::mul(Engine::mul(Engine::make_complex(Engine::mul(cos_Phi, Pi), zero), Dn), rn), Rho);

    // Ne1n
    // Ne1n[0] = cos_Phi * (n * n + n) * sin_Theta * Pi * rn / Rho / Rho;
    factor = Engine::mul(cos_Phi, Engine::mul(n_sq_plus_n, Engine::mul(sin_Theta, Pi)));
    Ne1n[0] = Engine::div(Engine::div(Engine::mul(Engine::make_complex(factor, zero), rn), Rho), Rho);

    // Ne1n[1] = cos_Phi * Tau * Dn * rn / Rho;
    Ne1n[1] = Engine::div(Engine::mul(Engine::mul(Engine::make_complex(Engine::mul(cos_Phi, Tau), zero), Dn), rn), Rho);

    // Ne1n[2] = -sin_Phi * Pi * Dn * rn / Rho;
    Ne1n[2] = Engine::div(Engine::mul(Engine::mul(Engine::make_complex(Engine::mul(Engine::sub(zero, sin_Phi), Pi), zero), Dn), rn), Rho);
  }

  static inline void calcPiTau(const RealV costheta, std::vector<RealV>& Pi, std::vector<RealV>& Tau) {
    int nmax = Pi.size();
    auto one = Engine::set(1.0);
    auto two = Engine::set(2.0);
    auto three = Engine::set(3.0);
    
    Pi[0] = one;
    Tau[0] = costheta;
    
    if (nmax > 1) {
      // Pi[1] = 3*costheta*Pi[0];
      Pi[1] = Engine::mul(three, Engine::mul(costheta, Pi[0]));
      // Tau[1] = 2*costheta*Pi[1] - 3*Pi[0];
      Tau[1] = Engine::sub(Engine::mul(two, Engine::mul(costheta, Pi[1])), Engine::mul(three, Pi[0]));
      
      for (int i = 2; i < nmax; i++) {
        auto i_val = Engine::set((FloatType)i);
        auto i_plus_1 = Engine::set((FloatType)(i + 1));
        auto i_plus_2 = Engine::set((FloatType)(i + 2));
        auto two_i_plus_1 = Engine::set((FloatType)(2 * i + 1));
        
        // Pi[i] = ((2*i + 1)*costheta*Pi[i - 1] - (i + 1)*Pi[i - 2])/i;
        auto term1 = Engine::mul(two_i_plus_1, Engine::mul(costheta, Pi[i - 1]));
        auto term2 = Engine::mul(i_plus_1, Pi[i - 2]);
        Pi[i] = Engine::div(Engine::sub(term1, term2), i_val);
        
        // Tau[i] = (i + 1)*costheta*Pi[i] - (i + 2)*Pi[i - 1];
        term1 = Engine::mul(i_plus_1, Engine::mul(costheta, Pi[i]));
        term2 = Engine::mul(i_plus_2, Pi[i - 1]);
        Tau[i] = Engine::sub(term1, term2);
      }
    }
  }
};

template <typename FloatType, typename Engine>
void RunFieldKernel(
    const std::vector<FloatType>& Xp,
    const std::vector<FloatType>& Yp,
    const std::vector<FloatType>& Zp,
    size_t start_index,
    size_t count,
    int nmax,
    const std::vector<FloatType>& size_param,
    const std::vector<std::complex<FloatType>>& refractive_index,
    const std::vector<std::vector<std::complex<FloatType>>>& aln,
    const std::vector<std::vector<std::complex<FloatType>>>& bln,
    const std::vector<std::vector<std::complex<FloatType>>>& cln,
    const std::vector<std::vector<std::complex<FloatType>>>& dln,
    MieBuffers<FloatType, Engine>& buffers,
    std::vector<std::vector<std::complex<FloatType>>>& Es,
    std::vector<std::vector<std::complex<FloatType>>>& Hs,
    std::vector<std::vector<FloatType>>& coords_polar
) {
    using RealV = typename Engine::RealV;
    using ComplexV = typename Engine::ComplexV;
    const size_t lanes = Engine::Lanes();

    // Buffers
    std::vector<FloatType> rho_buf(lanes), theta_buf(lanes), phi_buf(lanes);
    std::vector<FloatType> ml_re_buf(lanes), ml_im_buf(lanes);
    std::vector<int> layer_buf(lanes);
    std::vector<size_t> original_indices(lanes);

    // Loop over batches
    for (size_t i = 0; i < count; i += lanes) {
        size_t n_process = std::min(lanes, count - i);
        
        for (size_t j = 0; j < n_process; ++j) {
            size_t idx = start_index + i + j;
            FloatType x = Xp[idx];
            FloatType y = Yp[idx];
            FloatType z = Zp[idx];
            
            using nmm::sqrt; using nmm::acos; using nmm::atan2;
            FloatType rho = sqrt(x*x + y*y + z*z);
            FloatType theta = (rho > 0.0) ? acos(z/rho) : 0.0;
            FloatType phi = atan2(y, x);
            if (rho < 1e-5) rho = 1e-5;
            
            rho_buf[j] = rho;
            theta_buf[j] = theta;
            phi_buf[j] = phi;
            
            coords_polar[idx] = {rho, theta, phi};
            
            // GetIndexAtRadius logic
            unsigned int l = 0;
            std::complex<FloatType> ml;
            if (rho > size_param.back()) {
                l = size_param.size();
                ml = std::complex<FloatType>(1.0, 0.0);
            } else {
                for (int k = size_param.size() - 1; k >= 0 ; k--) {
                    if (rho <= size_param[k]) {
                        l = k;
                    }
                }
                ml = refractive_index[l];
            }
            
            ml_re_buf[j] = ml.real();
            ml_im_buf[j] = ml.imag();
            layer_buf[j] = l;
            original_indices[j] = idx;
        }
        
        // Fill remaining lanes with dummy data (copy of first)
        for (size_t j = n_process; j < lanes; ++j) {
            rho_buf[j] = rho_buf[0];
            theta_buf[j] = theta_buf[0];
            phi_buf[j] = phi_buf[0];
            ml_re_buf[j] = ml_re_buf[0];
            ml_im_buf[j] = ml_im_buf[0];
            layer_buf[j] = layer_buf[0];
        }

        RealV rho_v = Engine::load(rho_buf.data());
        RealV theta_v = Engine::load(theta_buf.data());
        RealV phi_v = Engine::load(phi_buf.data());
        ComplexV ml_v = Engine::make_complex(Engine::load(ml_re_buf.data()), Engine::load(ml_im_buf.data()));
        
        ComplexV rho_c = Engine::make_complex(rho_v, Engine::set(0.0));
        ComplexV z = Engine::mul(rho_c, ml_v);
        
        auto& D1 = buffers.D1;
        auto& D3 = buffers.D3;
        auto& Psi = buffers.Psi;
        auto& Zeta = buffers.Zeta;
        auto& PsiZeta = buffers.PsiZeta;
        
        evalDownwardD1<FloatType, Engine>(z, D1);
        evalUpwardPsi<FloatType, Engine>(z, D1, Psi);
        evalUpwardD3<FloatType, Engine>(z, D1, D3, PsiZeta);
        
        for(int k=0; k<=nmax; ++k) {
           auto psi_zeta_k = Engine::load(&PsiZeta[k * lanes]);
           auto psi_k = Engine::load(&Psi[k * lanes]);
           auto zeta_k = Engine::div(psi_zeta_k, psi_k);
           Engine::store(zeta_k, &Zeta[k * lanes]);
        }
        
        std::vector<RealV> Pi(nmax), Tau(nmax);
        NearFieldKernelHelper<FloatType, Engine>::calcPiTau(Engine::cos(theta_v), Pi, Tau);
        
        ComplexV E_v[3], H_v[3];
        auto zero = Engine::set(0.0);
        auto c_zero = Engine::make_complex(zero, zero);
        for(int k=0; k<3; ++k) { E_v[k] = c_zero; H_v[k] = c_zero; }
        
        auto c_i = Engine::make_complex(zero, Engine::set(1.0));
        auto c_one = Engine::make_complex(Engine::set(1.0), zero);
        
        // ipow
        ComplexV ipow[4] = {c_one, c_i, Engine::make_complex(Engine::set(-1.0), zero), Engine::make_complex(zero, Engine::set(-1.0))};

        // Pre-gather coefficients to improve cache locality
        // We transpose from [layer][n] to [n][lane]
        std::vector<FloatType> aln_re_all(nmax * lanes),
            aln_im_all(nmax * lanes);
        std::vector<FloatType> bln_re_all(nmax * lanes),
            bln_im_all(nmax * lanes);
        std::vector<FloatType> cln_re_all(nmax * lanes),
            cln_im_all(nmax * lanes);
        std::vector<FloatType> dln_re_all(nmax * lanes),
            dln_im_all(nmax * lanes);

        for (size_t j = 0; j < lanes; ++j) {
          int l = layer_buf[j];
          // Optimization: if multiple lanes have same layer, we could copy?
          // But simple loop is fine for now.
          for (int n = 0; n < nmax; ++n) {
            size_t idx = n * lanes + j;
            aln_re_all[idx] = aln[l][n].real();
            aln_im_all[idx] = aln[l][n].imag();
            bln_re_all[idx] = bln[l][n].real();
            bln_im_all[idx] = bln[l][n].imag();
            cln_re_all[idx] = cln[l][n].real();
            cln_im_all[idx] = cln[l][n].imag();
            dln_re_all[idx] = dln[l][n].real();
            dln_im_all[idx] = dln[l][n].imag();
          }
        }

        for (int n = 0; n < nmax; ++n) {
          int n1 = n + 1;
          auto rn = Engine::make_complex(Engine::set((FloatType)n1), zero);

          ComplexV M1o1n[3], M1e1n[3], N1o1n[3], N1e1n[3];
          ComplexV M3o1n[3], M3e1n[3], N3o1n[3], N3e1n[3];

          NearFieldKernelHelper<FloatType, Engine>::calcSpherHarm(
              z, theta_v, phi_v, Engine::load(&Psi[n1 * lanes]),
              Engine::load(&D1[n1 * lanes]), Pi[n], Tau[n],
              Engine::get_real(rn), M1o1n, M1e1n, N1o1n, N1e1n);
          NearFieldKernelHelper<FloatType, Engine>::calcSpherHarm(
              z, theta_v, phi_v, Engine::load(&Zeta[n1 * lanes]),
              Engine::load(&D3[n1 * lanes]), Pi[n], Tau[n],
              Engine::get_real(rn), M3o1n, M3e1n, N3o1n, N3e1n);

          auto En = Engine::mul(
              ipow[n1 % 4],
              Engine::make_complex(
                  Engine::div(Engine::set((FloatType)(2 * n1 + 1)),
                              Engine::set((FloatType)(n1 * n1 + n1))),
                  zero));

          // Load pre-gathered coefficients
          ComplexV aln_v =
              Engine::make_complex(Engine::load(&aln_re_all[n * lanes]),
                                   Engine::load(&aln_im_all[n * lanes]));
          ComplexV bln_v =
              Engine::make_complex(Engine::load(&bln_re_all[n * lanes]),
                                   Engine::load(&bln_im_all[n * lanes]));
          ComplexV cln_v =
              Engine::make_complex(Engine::load(&cln_re_all[n * lanes]),
                                   Engine::load(&cln_im_all[n * lanes]));
          ComplexV dln_v =
              Engine::make_complex(Engine::load(&dln_re_all[n * lanes]),
                                   Engine::load(&dln_im_all[n * lanes]));

          for(int k=0; k<3; ++k) {
            // Ediff = En*(cln*M1o1n[i] - c_i*dln*N1e1n[i] + c_i*aln*N3e1n[i] - bln*M3o1n[i]);
            auto term1 = Engine::mul(cln_v, M1o1n[k]);
            auto term2 = Engine::mul(c_i, Engine::mul(dln_v, N1e1n[k]));
            auto term3 = Engine::mul(c_i, Engine::mul(aln_v, N3e1n[k]));
            auto term4 = Engine::mul(bln_v, M3o1n[k]);
            auto Ediff = Engine::mul(En, Engine::sub(Engine::add(Engine::sub(term1, term2), term3), term4));
            E_v[k] = Engine::add(E_v[k], Ediff);
            
            // Hdiff = En*(-dln*M1e1n[i] - c_i*cln*N1o1n[i] + c_i*bln*N3o1n[i] + aln*M3e1n[i]);
            term1 = Engine::mul(Engine::sub(c_zero, dln_v), M1e1n[k]);
            term2 = Engine::mul(c_i, Engine::mul(cln_v, N1o1n[k]));
            term3 = Engine::mul(c_i, Engine::mul(bln_v, N3o1n[k]));
            term4 = Engine::mul(aln_v, M3e1n[k]);
            auto Hdiff = Engine::mul(En, Engine::add(Engine::add(Engine::sub(term1, term2), term3), term4));
            H_v[k] = Engine::add(H_v[k], Hdiff);
          }
        }

        // Add incident field (Vectorized)
        std::vector<FloatType> layer_f(lanes);
        for(size_t j=0; j<lanes; ++j) layer_f[j] = static_cast<FloatType>(layer_buf[j]);
        auto layer_v = Engine::load(layer_f.data());
        auto outer_layer_val = Engine::set(static_cast<FloatType>(refractive_index.size()));
        auto mask_incident = Engine::eq(layer_v, outer_layer_val);
        
        auto z_val = Engine::mul(rho_v, Engine::cos(theta_v));
        auto Ex = Engine::make_complex(Engine::cos(z_val), Engine::sin(z_val));
        
        auto sin_theta = Engine::sin(theta_v);
        auto cos_theta = Engine::cos(theta_v);
        auto sin_phi = Engine::sin(phi_v);
        auto cos_phi = Engine::cos(phi_v);
        
        auto E_inc_rho = Engine::mul(Ex, Engine::mul(cos_phi, sin_theta));
        auto E_inc_theta = Engine::mul(Ex, Engine::mul(cos_phi, cos_theta));
        auto E_inc_phi = Engine::mul(Ex, Engine::sub(zero, sin_phi));
        
        E_v[0] = Engine::add(E_v[0], Engine::select(mask_incident, E_inc_rho, c_zero));
        E_v[1] = Engine::add(E_v[1], Engine::select(mask_incident, E_inc_theta, c_zero));
        E_v[2] = Engine::add(E_v[2], Engine::select(mask_incident, E_inc_phi, c_zero));
        
        auto Hy = Ex;
        auto H_inc_rho = Engine::mul(Hy, Engine::mul(sin_theta, sin_phi));
        auto H_inc_theta = Engine::mul(Hy, Engine::mul(cos_theta, sin_phi));
        auto H_inc_phi = Engine::mul(Hy, cos_phi);
        
        H_v[0] = Engine::add(H_v[0], Engine::select(mask_incident, H_inc_rho, c_zero));
        H_v[1] = Engine::add(H_v[1], Engine::select(mask_incident, H_inc_theta, c_zero));
        H_v[2] = Engine::add(H_v[2], Engine::select(mask_incident, H_inc_phi, c_zero));

        // Apply magnetic field factor
        auto hffact_v = Engine::div(ml_v, Engine::set(static_cast<FloatType>(nmie::cc_ * nmie::mu_)));
        for(int k=0; k<3; ++k) {
            H_v[k] = Engine::mul(H_v[k], hffact_v);
        }

        // Store results
        std::vector<FloatType> E_re(lanes), E_im(lanes);
        std::vector<FloatType> H_re(lanes), H_im(lanes);
        
        for(int k=0; k<3; ++k) {
          Engine::store(Engine::get_real(E_v[k]), E_re.data());
          Engine::store(Engine::get_imag(E_v[k]), E_im.data());
          Engine::store(Engine::get_real(H_v[k]), H_re.data());
          Engine::store(Engine::get_imag(H_v[k]), H_im.data());
          
          for(size_t j=0; j<n_process; ++j) {
            Es[original_indices[j]][k] = std::complex<FloatType>(E_re[j], E_im[j]);
            Hs[original_indices[j]][k] = std::complex<FloatType>(H_re[j], H_im[j]);
          }
        }
    }
}
}

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
    aln_.clear(); bln_.clear(); cln_.clear(); dln_.clear();

    std::complex<FloatType> c_one(1.0, 0.0), c_zero(0.0, 0.0);

    const int L = refractive_index_.size();

    aln_.resize(L + 1);
    bln_.resize(L + 1);
    cln_.resize(L + 1);
    dln_.resize(L + 1);
    for (int l = 0; l <= L; l++) {
      aln_[l].resize(nmax_, static_cast<FloatType>(0.0));
      bln_[l].resize(nmax_, static_cast<FloatType>(0.0));
      cln_[l].resize(nmax_, static_cast<FloatType>(0.0));
      dln_[l].resize(nmax_, static_cast<FloatType>(0.0));
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

    int print_precision = 16;
#ifdef MULTI_PRECISION
    print_precision = MULTI_PRECISION;
#endif
    // Check the result and change  aln_[0][n] and aln_[0][n] for exact zero
    int print_count = 0;
    for (int n = 0; n < nmax_; ++n) {
      if (cabs(aln_[0][n]) > 1e-10 && print_count < 2)  {
        print_count++;
        std::cout<< std::setprecision(print_precision)
                 << "Warning: Potentially unstable calculation of aln[0]["
                 << n << "] = "<< aln_[0][n] << " which is expected to be exact zero!"<<std::endl;
      }
      if (cabs(bln_[0][n]) > 1e-10  && print_count < 2)  {
        print_count++;
        std::cout<< std::setprecision(print_precision)
                 << "Warning: Potentially unstable calculation of bln[0]["
                 << n << "] = "<< bln_[0][n] << " which is expected to be exact zero!" <<std::endl;
      }
      aln_[0][n] = 0.0;
      bln_[0][n] = 0.0;
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
  // Main trouble of near-field evaluations is supposed to originate from special functions
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
                                  std::vector<std::complex<evalType> > &H,
                                  std::vector<bool> &isConvergedE,
                                  std::vector<bool> &isConvergedH,
                                  bool isMarkUnconverged)  {
    auto nmax = Psi.size() - 1;
    std::complex<evalType> c_zero(0.0, 0.0), c_i(0.0, 1.0), c_one(1.0, 0.0);
//    auto c_nan = ConvertComplex<FloatType>(std::complex<double>(std::nan(""), std::nan("")));
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

    const unsigned L = refractive_index_.size();
    for (int n = 0; n < nmax_; n++) {
      cln_[L][n] = c_zero;
      dln_[L][n] = c_zero;
    }

    unsigned int l;
    GetIndexAtRadius(Rho, ml, l);

    isConvergedE = {false, false, false}, isConvergedH = {false, false, false};
//    evalType E0 = 0, H0=0;
    std::vector< std::complex<evalType> > Ediff_prev = {{0.,0.},{0.,0.},{0.,0.}},
        Hdiff_prev = {{0.,0.},{0.,0.},{0.,0.}};
    for (unsigned int n = 0; n < nmax; n++) {
      if ( isConvergedE[0] && isConvergedE[1] && isConvergedE[2]
          && isConvergedH[0] && isConvergedH[1] && isConvergedH[2]) {
        std::cout<<"Near-field early convergence at nmax = "<<n+1<<std::endl;
        break;
      }
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
      auto aln = ConvertComplex<evalType>(aln_[l][n]);
      auto bln = ConvertComplex<evalType>(bln_[l][n]);
      auto cln = ConvertComplex<evalType>(cln_[l][n]);
      auto dln = ConvertComplex<evalType>(dln_[l][n]);
      for (int i = 0; i < 3; i++) {
        if (isConvergedE[i] && isConvergedH[i]) continue; // TODO is it safe?
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
        if (n>0) {
          if (
              (cabs(Ediff_prev[i]) <= cabs(E[i]) * nearfield_convergence_threshold_)
                  &&  (cabs(Ediff) <= cabs(E[i]) * nearfield_convergence_threshold_)
              ) isConvergedE[i] = true;
          if (
              (cabs(Hdiff_prev[i]) <= cabs(H[i]) * nearfield_convergence_threshold_)
                  &&  (cabs(Hdiff) <= cabs(H[i]) * nearfield_convergence_threshold_)
              ) isConvergedH[i] = true;
        }
        Ediff_prev[i] = Ediff;
        Hdiff_prev[i] = Hdiff;

        if ((!isConvergedH[i] || !isConvergedE[i]) && n==nmax-1 && GetFieldConvergence()) {
          std::cout<<"Econv:"<<cabs(Ediff)/cabs(E[i])<<" Hconv:"<<cabs(Hdiff)/cabs(H[i])<<std::endl;

        }
        if (mode_n_ == Modes::kAll) {
          // electric field E [V m - 1] = EF*E0
          E[i] += Ediff;
          H[i] += Hdiff;
          continue;
        }
        if (n == 0) {

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

    // Add the incident field
    if(l==L) {
      const auto z = Rho*cos_t(Theta);
      const auto Ex = std::complex<evalType>(cos_t(z), sin_t(z));
      E[0] +=  Ex*cos_t(Phi)*sin_t(Theta);
      E[1] +=  Ex*cos_t(Phi)*cos_t(Theta);
      E[2] += -Ex*sin_t(Phi);
      const auto Hy = Ex;
      H[0] += Hy*sin_t(Theta)*sin_t(Phi);
      H[1] += Hy*cos_t(Theta)*sin_t(Phi);
      H[2] += Hy*cos_t(Phi);
    }

    if( (!isConvergedE[0] || !isConvergedE[1] ||!isConvergedE[2] ||
        !isConvergedH[0] || !isConvergedH[1] ||!isConvergedH[2] ) && GetFieldConvergence()) {
      std::cout << "Field evaluation failed to converge an nmax = "<< nmax << std::endl;
      std::cout << "Near-field convergence threshold: "<<nearfield_convergence_threshold_<<std::endl;
      if (isMarkUnconverged) {  //mark as NaN
        for(auto &ee :E) ee /= c_zero;
        for(auto &ee :H) ee /= c_zero;
      }
    }

    // magnetic field
    std::complex<evalType> hffact = ml/static_cast<evalType>(nmie::cc_*nmie::mu_);
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
  void MultiLayerMie<FloatType>::RunFieldCalculation(bool isMarkUnconverged) {
    (void)isMarkUnconverged;
    // Calculate scattering coefficients an_ and bn_
    calcScattCoeffs();
    // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
    calcExpanCoeffs();

    // Zero out cln and dln for the outer layer (L) to avoid adding incident field twice
    const unsigned L = refractive_index_.size();
    std::complex<FloatType> c_zero_host(0.0, 0.0);
    for (int n = 0; n < nmax_; n++) {
      cln_[L][n] = c_zero_host;
      dln_[L][n] = c_zero_host;
    }

    isConvergedE_ = {true, true, true}, isConvergedH_ = {true, true, true};
    Es_.clear(); Hs_.clear(); coords_polar_.clear();
    long total_points = coords_[0].size();
    
    Es_.resize(total_points);
    Hs_.resize(total_points);
    coords_polar_.resize(total_points);
    for(long i=0; i<total_points; ++i) {
        Es_[i].resize(3);
        Hs_[i].resize(3);
    }

#ifdef WITH_HWY
    using Engine = HighwayEngine<FloatType>;
#else
    using Engine = ScalarEngine<FloatType>;
#endif

    MieBuffers<FloatType, Engine> buffers;
    buffers.resize(nmax_, 1);

    RunFieldKernel<FloatType, Engine>(
        coords_[0], coords_[1], coords_[2],
        0, total_points,
        nmax_,
        size_param_,
        refractive_index_,
        aln_, bln_, cln_, dln_,
        buffers,
        Es_, Hs_,
        coords_polar_
    );

    convertFieldsFromSphericalToCartesian();
  }  //  end of MultiLayerMie::RunFieldCalculation()

// TODO do we really need this eval_delta()?
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
void MultiLayerMie<FloatType>::calcMieSeriesNeededToConverge(const FloatType Rho, int nmax_in) {
  if (nmax_in < 1) {
    auto required_near_field_nmax = calcNmax(Rho);
    SetMaxTerms(required_near_field_nmax);
  } else {
    SetMaxTerms(nmax_in);
  }
  // Calculate scattering coefficients an_ and bn_
  calcScattCoeffs();
  // We might be limited with available machine precision
  available_maximal_nmax_ = nmax_;
  // Calculate expansion coefficients aln_,  bln_, cln_, and dln_
  calcExpanCoeffs();
}


template <typename FloatType>
void MultiLayerMie<FloatType>::calcRadialOnlyDependantFunctions(const double from_Rho, const double to_Rho,
                                                                std::vector<std::vector<std::complex<FloatType> > > &Psi,
                                                                std::vector<std::vector<std::complex<FloatType> > > &D1n,
                                                                std::vector<std::vector<std::complex<FloatType> > > &Zeta,
                                                                std::vector<std::vector<std::complex<FloatType> > > &D3n,
                                                                int nmax_in) {
  auto radius_points = Psi.size();
  std::vector<std::vector<std::complex<FloatType> > > PsiZeta(radius_points);
  double delta_Rho = eval_delta<double>(radius_points, from_Rho, to_Rho);
  for (unsigned int j=0; j < radius_points; j++) {
    auto Rho = static_cast<FloatType>(from_Rho + j*delta_Rho);
//    if (Rho < 1e-5) Rho = 1e-5; // TODO do we need this?.
    int near_field_nmax = nmax_in;
    if (nmax_in < 1) near_field_nmax = calcNmax(Rho);

    // Skip if not enough terms in Mie series (i.e. required near field nmax > available terms )
    if (near_field_nmax > available_maximal_nmax_)  near_field_nmax = available_maximal_nmax_;
    Psi[j].resize(near_field_nmax + 1, static_cast<FloatType>(0.0)); D1n[j].resize(near_field_nmax + 1, static_cast<FloatType>(0.0));
    Zeta[j].resize(near_field_nmax + 1, static_cast<FloatType>(0.0)); D3n[j].resize(near_field_nmax + 1, static_cast<FloatType>(0.0));
    PsiZeta[j].resize(near_field_nmax + 1, static_cast<FloatType>(0.0));
    std::complex<FloatType> ml;
    GetIndexAtRadius(Rho, ml);
    auto z = Rho*ml;
    evalDownwardD1<FloatType>(z, D1n[j]);
    evalUpwardPsi<FloatType>(z,  D1n[j], Psi[j]);
    evalUpwardD3<FloatType> (z, D1n[j], D3n[j], PsiZeta[j]);
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
                                                        const bool isMarkUnconverged,
                                                        int nmax_in) {
  if (from_Rho > to_Rho || from_Theta > to_Theta || from_Phi > to_Phi
      || outer_arc_points < 1 || radius_points < 1
      || from_Rho < 0.)
    throw std::invalid_argument("Error! Invalid argument for RunFieldCalculationPolar() !");
//  auto nmax_old = nmax_;
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
  calcMieSeriesNeededToConverge(to_Rho, nmax_in);

  std::vector<std::vector<FloatType> >  Pi(theta_points), Tau(theta_points);
  calcPiTauAllTheta(from_Theta, to_Theta, Pi, Tau);

  std::vector<std::vector<std::complex<FloatType> > > Psi(radius_points), D1n(radius_points),
      Zeta(radius_points), D3n(radius_points), PsiZeta(radius_points);
  calcRadialOnlyDependantFunctions(from_Rho, to_Rho,
                                   Psi, D1n, Zeta, D3n,
                                   nmax_in);

//  std::cout<<"Done evaluation of special functions."<<std::endl;
  double delta_Rho = eval_delta<double>(radius_points, from_Rho, to_Rho);
  double delta_Theta = eval_delta<double>(theta_points, from_Theta, to_Theta);
  double delta_Phi = eval_delta<double>(phi_points, from_Phi, to_Phi);
  Es_.clear(); Hs_.clear(); coords_polar_.clear();
  std::vector<bool> isConvergedE = {false, false, false}, isConvergedH = {false, false, false};
  isConvergedE_ = {true, true, true}, isConvergedH_ = {true, true, true};
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
        std::vector<std::complex<double> > Es(3), Hs(3);
        calcFieldByComponents( Rho, Theta, Phi,
                               Psi_dp, D1n_dp, Zeta_dp, D3n_dp, Pi_dp, Tau_dp,
                               Es, Hs, isConvergedE, isConvergedH,
                               isMarkUnconverged);
        UpdateConvergenceStatus(isConvergedE, isConvergedH);
        Es_.push_back(ConvertComplexVector<FloatType>(Es));
        Hs_.push_back(ConvertComplexVector<FloatType>(Hs));
      }
    }
  }
  convertFieldsFromSphericalToCartesian();
}


template <typename FloatType>
void MultiLayerMie<FloatType>::UpdateConvergenceStatus(std::vector<bool> isConvergedE, std::vector<bool> isConvergedH) {
  for (int i = 0; i< 3; i++) isConvergedE_[i] = isConvergedE_[i] && isConvergedE[i];
  for (int i = 0; i< 3; i++) isConvergedH_[i] = isConvergedH_[i] && isConvergedH[i];
}


template <typename FloatType>
bool MultiLayerMie<FloatType>::GetFieldConvergence () {
  bool convergence = true;
  for (auto conv:isConvergedE_) convergence = convergence && conv;
  for (auto conv:isConvergedH_) convergence = convergence && conv;
  return convergence;
}

template <typename FloatType>
void MultiLayerMie<FloatType>::RunFieldCalculationCartesian(const int first_side_points,
                                                            const int second_side_points,
                                                            const double relative_side_length,
                                                            const int plane_selected,
                                                            const double at_x, const double at_y,
                                                            const double at_z,
                                                            const bool isMarkUnconverged,
                                                            const int nmax_in) {
  SetMaxTerms(nmax_in);
  std::vector<FloatType> Xp(0), Yp(0), Zp(0);
  if (size_param_.size()<1) throw "Expect size_param_ to have at least one element before running a simulation";
  const FloatType total_R = size_param_.back();
  const FloatType second_side_max_coord_value = total_R * relative_side_length;
  // TODO add test if side_1_points <= 1 or side_2_points <= 1
  const FloatType space_step = second_side_max_coord_value*2/( (second_side_points<2 ? 2 : second_side_points) - 1.0);
  auto push_coords = [&](const int nx, const int ny, const int nz) {
    const FloatType xi = at_x*total_R - space_step*(nx-1)/2;
    const FloatType yi = at_y*total_R - space_step*(ny-1)/2;
    const FloatType zi = at_z*total_R - space_step*(nz-1)/2;
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
          Xp.push_back(xi + static_cast<FloatType>(i) * space_step);
          Yp.push_back(yi + static_cast<FloatType>(j) * space_step);
          Zp.push_back(zi + static_cast<FloatType>(k) * space_step);
        }
      }
    }
  };
  // TODO add test to check that side_2_points is for z-axis
  if (plane_selected == Planes::kEk) push_coords(first_side_points, 1, second_side_points);
  if (plane_selected == Planes::kHk) push_coords(1, first_side_points, second_side_points);
  if (plane_selected == Planes::kEH) push_coords(first_side_points, second_side_points, 1);
  const unsigned int total_size = first_side_points*second_side_points;
  if (Xp.size() != total_size || Yp.size() != total_size || Zp.size() != total_size)
    throw std::invalid_argument("Error! Wrong dimension of field monitor points for cartesian grid!");
  SetFieldCoords({Xp, Yp, Zp});
  RunFieldCalculation(isMarkUnconverged);
}  // end of void MultiLayerMie<FloatType>::RunFieldCalculationCartesian(...)

}  // end of namespace nmie
#endif  // SRC_NMIE_NEARFIELD_HPP_
