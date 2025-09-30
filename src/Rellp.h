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

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELLP_H_
#define RELLP_H_

#include "Relbase.h"

extern "C" {
#include "relutility.h"
#include "reltable.h"
#include "common.h"
}

#define NHBINS_VERTICALLY_EXTENDED_SOURCE 50

typedef struct {

  double *heightArr;  //nh+1 bins
  double *heightMean;
  double *beta;
  int nh;

} extPrimSource;


emisProfile *calc_emis_profile(double *r, int nr, const relParam *param, int *status);

void get_emis_bkn(double *emis, const double *re, int nr,
                  double index1, double index2, double rbr);

// relxill_emis - relxill3 and relxill_emis emissivity functions
void get_emis_bkn3(double *emis, const double *re, int nr,
                  double index1, double index2, double index3, double rbr, double rbr2);

void get_emis_rbin(double *emis, const double *re, int nr, int num_emis, double rin, double rout,
                  double index1, double index2, double index3, double index4, double index5,
                  double index6, double index7, double index8, double index9, double index10,
                  double index11, double index12, double index13, double index14, double index15,
                  double index16, double index17, double index18, double index19, double index20
                  );
// -- end

void calc_emis_jet_point_source(emisProfile *emis_profile, const relParam *param, double height, double beta,
                                       lpTable *tab, int *status);

lpTable* get_lp_table(int* status);

void get_emis_jet(emisProfile *, const relParam *param, int *status);

void rebin_emisprofile_on_radial_grid(emisProfile *emis_prof, const emisProfile* emis_prof_tab, int *status);

void apply_emis_fluxboost_source_disk(emisProfile *emisProf, double a, double height, double gamma, double beta);

int modelLampPostPointsource(const relParam *param);

extPrimSource *getExtendedJetGeom(const relParam *param, int *status);

////////////

void free_cached_lpTable(void);

lpReflFrac *new_lpReflFrac(int *status);
void free_lpReflFrac(lpReflFrac **str);

emisProfile *new_emisProfile(double *re, int nr, int *status);
void free_emisProfile(emisProfile *emis_profile);

extPrimSource *new_extendedPrimarySource(int nh, int *status);
void free_extendedPrimarySource(extPrimSource *source);

int is_emis_grid_ascending(const emisProfile* emis);

#endif /* RELLP_H_ */
