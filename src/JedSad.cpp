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

#include "JedSad.h"

extern "C" {
#include "xilltable.h"
#include "relutility.h"
}

#define CHECK_STATUS_CFITSIO(status) \
 if (EXIT_SUCCESS!=status) {printf(" error in cfitsio function with code %i\n",status); throw std::exception();}



void JedsadTable::read_jstable_params(fitsfile* fptr){

  int extver = 0;
  string extname = "PARAMETERS";
  int status = EXIT_SUCCESS;
  fits_movenam_hdu_cpp(fptr, BINARY_TBL, extname.c_str(), extver, &status);
  CHECK_STATUS_CFITSIO(status)

  // we know the column numbers
  int colnum_n = 2;
  int colnum_vals = 3;

  long n;
  fits_get_num_rows(fptr, &n, &status);
  CHECK_STATUS_CFITSIO(status)

  if (num_params() != n) {
    printf("wrong format of the JEDSAD table (not the correct number of parameters)");
    throw std::exception();
  }

  int anynul = 0;
  double nullval = 0.0;

  int ii;
  char strnull[10];
  strcpy(strnull, " ");

  fits_read_col(fptr, TINT, colnum_n, 1, 1, n, &nullval, num_param_vals, &anynul, &status);
  CHECK_STATUS_CFITSIO(status)

  assert(param_vals.size() == n);
  for (ii = 0; ii < n; ii++) {
    /** the we load the parameter values **/
    param_vals[ii] = new double[num_param_vals[ii]];
    fits_read_col(fptr, TDOUBLE, colnum_vals, ii + 1, 1, num_param_vals[ii], &nullval, param_vals[ii], &anynul,&status);
    CHECK_STATUS_CFITSIO(status)
  }

  int n_data_entries = 1;
  for (ii = 0; ii < n; ii++){
    n_data_entries *= num_param_vals[ii];
  }

  // read the data
  string extname_data = "DATA";
  fits_movenam_hdu_cpp(fptr, BINARY_TBL, extname_data.c_str(), extver, &status);
  CHECK_STATUS_CFITSIO(status)

  int n_columns_data;
  fits_get_num_cols(fptr, &n_columns_data, &status);
  CHECK_STATUS_CFITSIO(status)
  assert(n_columns_data == data.size());

  long n_rows_data;
  fits_get_num_rows(fptr, &n_rows_data, &status);
  CHECK_STATUS_CFITSIO(status)
  assert(n_rows_data == n_data_entries);

  for(ii = 0; ii< n_columns_data; ii++) {
    data[ii] = new double[n_rows_data];
    fits_read_col(fptr, TDOUBLE, ii+1, 1, 1, n_rows_data, &nullval, data[ii], &anynul, &status);
  }

}



void JedsadTable::read_table() {


  // init sizes
  param_vals.resize(num_params());
  data.resize(JedsadTableInformation::instance().n_data);

  int status = EXIT_SUCCESS;
  fitsfile *fptr = open_fits_table_stdpath(m_fullfilename.c_str(), &status);
  CHECK_STATUS_CFITSIO(status);

  assert(fptr != nullptr);
  read_jstable_params(fptr);

}

void JedsadTable::interpolate(const double *param_array) {

  if (data.empty()) {
    read_table();
  }

}

PrimespecParams convert_jedsad_to_primaryspec_params(const ModelDefinition& model_definition){

  auto jedsadParams = JedsadParams(model_definition);

  double* param_array = jedsadParams.get_param_array();
  JedsadTable::instance().interpolate(param_array);
  free(param_array);



  return PrimespecParams(-1.0, -1.0);
}


 // jedsadTable *cached_jedsad_table = nullptr;

