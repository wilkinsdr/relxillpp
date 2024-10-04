/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2024 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef JEDSAD_H_
#define JEDSAD_H_

#include <utility>
#include <array>
#include <numeric>
#include "ModelDefinition.h"

extern "C" {
#include "relutility.h"
}

enum JSPar{
  mass,
  rj,
  mdot,
  risco,
  omega,
  mu,
  ms,
  b,
  p
};


class JedsadTableInformation{

 public:
  /**
 * use this instance to always and securely access the database, without
 * needing to initialize it anywhere else
 */
  static JedsadTableInformation &instance() {
    static auto *instance = new JedsadTableInformation();
    return *instance;
  }

  // this is the order and name of the parameters in the table
  const std::vector<JSPar> param_names_table = {
      JSPar::mass,
      JSPar::rj,
      JSPar::mdot,
      JSPar::risco,
      JSPar::omega,
      JSPar::mu,
      JSPar::ms,
      JSPar::b,
      JSPar::p
  };

  string get_filename(){
    return m_filename;
  }

  const int n_data = 2;
  const int n_params = 9;

 private:
  const string m_filename = "Refl_XrB_DauserTest_lb3.fits";  // TODO: make it readable by env variable

};


class JedsadParams{

 public:
  explicit JedsadParams(const ModelDefinition& _model_params){
    jedsad_params[JSPar::mass] = _model_params[XPar::mass];
    jedsad_params[JSPar::rj] = _model_params[XPar::rj];
    jedsad_params[JSPar::mdot] = _model_params[XPar::mdot];
    jedsad_params[JSPar::risco] = _model_params[XPar::rin]; // need to make sure it is converted to the correct unit
    jedsad_params[JSPar::omega] = _model_params[XPar::omega];
    jedsad_params[JSPar::mu] = _model_params[XPar::mu];
    jedsad_params[JSPar::ms] = _model_params[XPar::ms];
    jedsad_params[JSPar::b] = _model_params[XPar::b];
    jedsad_params[JSPar::p] = _model_params[XPar::p];
  }

  double* get_param_array(){
    auto array = new double[jedsad_params.size()];
    for (int ii=0; ii<jedsad_params.size(); ii++){
      // use the correct order as defined in the Table (in JedsadTableInformation)
      array[ii] = jedsad_params[JedsadTableInformation::instance().param_names_table[ii]];
    }
    return array;
  }

  std::unordered_map<JSPar, double> jedsad_params;
};



class JedsadTable{

 public:
  JedsadTable(){
    m_fullfilename = JedsadTableInformation::instance().get_filename();  // the standard filename
    read_table();
  };

  size_t num_params() {
    return param_names_table.size();
  };

  static JedsadTable * const instance;
  void read_table();
  void interpolate(const double* param_array);

  int *num_param_vals = new int[num_params()];
  std::vector<double*> param_vals;
  std::vector<double*> data;

  void update_filename(){
    // check env variable (TBD)
    // if different, update it
  }

 private:
  string m_fullfilename;
  const std::vector<JSPar> param_names_table = JedsadTableInformation::instance().param_names_table;
  void read_jstable_params(fitsfile *fptr);
} ;


class PrimespecParams{ ;

 public:
  PrimespecParams(double _gamma, double _ecut) : gamma{_gamma}, ecut{_ecut} { };

  double ecut;
  double gamma;
};


PrimespecParams convert_jedsad_to_primaryspec_params(const ModelDefinition& model_definition);


#endif /* JEDSAD_H_ */
