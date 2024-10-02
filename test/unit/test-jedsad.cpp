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

    Copyright 2021 Thomas Dauser, Remeis Observatory & ECAP
*/

#include "catch2/catch_amalgamated.hpp"
#include "LocalModel.h"
#include "xspec_wrapper_lmodels.h"
#include "XspecSpectrum.h"
#include "common-functions.h"


TEST_CASE(" Execute relxill jedsad", "[jedsad]") {
  DefaultSpec default_spec{};

  LocalModel lmod(ModelName::relxill_jedsad);

  auto spec = default_spec.get_xspec_spectrum();
  REQUIRE_NOTHROW(lmod.eval_model(spec));

  double sum = sum_flux(spec.flux, spec.num_flux_bins());
  REQUIRE(sum > 1e-8);

}

