#ifndef SRC_READ_SPECTRA_READ_SPECTRA_H_
#define SRC_READ_SPECTRA_READ_SPECTRA_H_
/**
 * @file   read-spectra.h
 * @author Konstantin Ladutenko <kostyfisik at gmail (.) com>
 * @date   Wed Mar 11 11:19:34 2015
 * @copyright 2015 Konstantin Ladutenko
 *
 * @brief  Read complex spectra from file in format 'WL real imag'
 * 
 * read-spectra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * read-spectra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with read-spectra.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */
#include <vector>
#include <string>
#include <complex>
#include <utility>
namespace read_spectra {
  class ReadSpectra {  // will throw for any error
   public:
    ReadSpectra& ReadFromFile(std::string filename);
    ReadSpectra& ResizeToComplex(double from_wl, double to_wl, int samples);
    ReadSpectra& ToIndex();
    std::complex<double> at(double wavelength);
    void PrintData();    
    std::vector< std::pair< double, std::complex<double> > >&
      GetIndex(){return data_complex_index_;};
  private:
    std::vector< std::vector<double> > data_;
    std::vector< std::pair< double, std::complex<double> > > data_complex_;
    std::vector< std::pair< double, std::complex<double> > > data_complex_index_;
    void PermittivityToIndex();
  };  // end of class ReadSpectra
}  // end of namespase read_spectra
#endif  // SRC_READ_SPECTRA_READ_SPECTRA_H_
