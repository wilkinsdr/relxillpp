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

#include "Relphysics.h"
#include "Rellp_Extended.h" // for tabulated energy shifts only

extern "C" {
#include "common.h"
#include "relutility.h"
}

// calculate the u^t component of the 4-velocity of a thin accretion disk
// (see Bardeen, Press, Teukolsky; 1972)
double ut_disk(double r, double a) {
  return ((r * sqrt(r) + a) /
      (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
}

double calc_proper_area_ring(double rlo, double rhi, double a) {

  double rmean = 0.5 * (rlo + rhi);

  assert(a < 1);
  assert(a >= -1);

  double rho2 = rmean * rmean + a * a;  // this is rho^2
  double del = rmean * rmean - 2 * rmean + a * a;

  double area_gr = 2 * M_PI / sqrt(del);
  //  area_gr *= sqrt((rho2 * rho2 + 2 * a * a * rmean));
  area_gr *= sqrt((rho2 * rho2 - a * a * del));

  area_gr *= (rhi - rlo);

  return area_gr;
}

// Black Body Spectrum (returns ph/s/cm²/keV)  (verified with Xspec bbody model)
// ToDO: use a properly integrated model (?)x
void bbody_spec(const double *ener, int n, double *spec, double temperature, double gfac) {
  for (int ii = 0; ii < n; ii++) {
    double emean = 0.5 * (ener[ii] + ener[ii + 1]);
    spec[ii] = pow((emean / gfac), 2) / (exp(emean / gfac / temperature) - 1);
  }
}

// calculate the temperature for a given radius r and respective Rin following SS73
static double disk_temperature_alpha(double r, double Rin) {
  return pow((r / Rin), (-3. / 4)) * pow((1 - sqrt(Rin / r)), (1. / 4));
}

// temperature of a disk as used by diskbb (equal to alpha disk for large radii)
static double disk_temperature_diskbb(double r, double Rin, double Tin) {
  return pow(Tin * (r / Rin), (-3. / 4));  // if we use Rin, we don't get T=Tin at the inner edge
}

// temperature profile for the given radii; use Tmax at Rmax following Poutanan+2007
void disk_Tprofile_alpha(const double *rlo, const double *rhi, double *temp, int n, double Rin, double Tin) {

  // use Tin as Tmax
  double Rmax = pow((1.5), (4. / 5)) * Rin;  // Poutanen+2007
  double K = Tin / disk_temperature_alpha(Rmax, Rin);

  double RminGrid = 0.5 * (rlo[0] + rhi[0]);
  assert(RminGrid >= Rin);

  for (int ii = 0; ii < n; ii++) {
    temp[ii] = K * disk_temperature_alpha(0.5 * (rlo[ii] + rhi[ii]), Rin);
  }
}

// temperature profile for the given radii; use Tmax at Rmax following Poutanan+2007
void disk_Tprofile_diskbb(const double *rlo, const double *rhi, double *temp, int n, double Tin) {

  double Rin = 0.5 * (rlo[0] + rhi[0]); // needs to be rmean, such that we really get Tin at the inner zone
  for (int ii = 0; ii < n; ii++) {
    temp[ii] = disk_temperature_diskbb(0.5 * (rlo[ii] + rhi[ii]), Rin, Tin);
  }
}

double *get_tprofile(double *rlo, double *rhi, const int nrad, double Rin, double Tin, int type, int *status) {

  CHECK_STATUS_RET(*status, nullptr);

  auto tprofile = new double[nrad];
  CHECK_MALLOC_RET_STATUS(tprofile, status, nullptr)

  if (type == TPROFILE_ALPHA) {
    disk_Tprofile_alpha(rlo, rhi, tprofile, nrad, Rin, Tin);
  } else if (type == TPROFILE_DISKBB) {
    disk_Tprofile_diskbb(rlo, rhi, tprofile, nrad, Tin);
  } else {
    RELXILL_ERROR(" failed getting the temperature profile ", status);
    printf("    reason: type=%i is unknown\n", type);
  }

  return tprofile;
}




/*** we calculate the disk density from  Shakura & Sunyaev (1973)
 *    - for zone A as describe in their publication,  formula (Eq 2.11)
 *    - only the radial dependence is picked up here  (viscosity alpha=const.)
 *    - normalized such that dens(rms) = 1
 *    - rms is given in units of Rg
 *                                                               ***/
double density_ss73_zone_a(double radius, double rms) {
  return pow( (radius/rms) , (3. / 2)) * pow((1 - sqrt(rms  / radius)), -2);
}

double relat_abberation(double del, double beta) {
  return acos((cos(del) - beta) / (1 - beta * cos(del)));
}

/**
 * @brief calculates the energy shift from a lamp post source at h to infinity
 */
double calc_g_inf(double height, double a) {
  return sqrt(1.0 - (2 * height / (height * height + a * a)));
}

/**
<<<<<<< HEAD
 * @brief calculates the energy shift from a ring source to infinity (NOT YET FULLY IMPLEMENTED)
=======
 * @brief calculates the time-dilation component of energy shift from a ring source to infinity
>>>>>>> extended_source_emis
 * @details cosinetheta is the solution of equation r_ring = sqrt(h^2 * tan(theta)^2 + a^2 * sin(theta)^2) for cos(theta)
 * Theta is polar angle of the ring, between 0 and pi/2 (r_ring^2 + height^2 are always > a^2)
 * @param height height of the ring source
 * @param r_ring radius of the ring
 * @param a spin of the black hole
 * @return energy shift (only e^nu component in this function)
 */
double calc_g_inf_ext(double height, double r_ring, double a) {
  double g_inf;
  if (r_ring <= 1.0e-6){ // lamp-post source
    g_inf = calc_g_inf(height, a);
  } else {
    if (height == 0.0) { // the disk case, cosinetheta = 0
      double r = sqrt(r_ring * r_ring + a * a);
      double Sigma = r * r;
      double Delta = r * r - 2.0 * r + a * a;
      double rho_squared = (r * r + a * a) * (r * r + a * a) - Delta * a * a;
      g_inf = sqrt(Sigma * Delta / rho_squared);
    } else {
      if (a == 0.0){ // special case of zero spin. The block below would give inf
        double cosine_theta = height / sqrt(height * height + r_ring * r_ring);
        double r = height / cosine_theta;
        g_inf = sqrt(1 - 2.0 / r);
      } else { // general case
        /// can be replaced with a function
        double cosine_theta = sqrt((sqrt(4.0 * a * a * height * height + (height * height +
            r_ring * r_ring - a * a) * (height * height + r_ring * r_ring - a * a)) +
                a * a - height * height - r_ring * r_ring) / (2.0 * a * a));
        double sine_theta = sqrt(1.0 - pow(cosine_theta, 2));
        double r = height / cosine_theta; // spherical radius of the ring, not cylindrical!
        double Sigma = r * r + a * a * cosine_theta * cosine_theta;
        double Delta = r * r - 2.0 * r + a * a;
        double rho_squared = (r * r + a * a) * (r * r + a * a) - Delta * a * a * sine_theta * sine_theta;
        //double omega = 2.0 * a * r / rho_squared;
        // previous version returned just the lamp-post shift, which led to stronger redshift of the ring primary source
        //return sqrt(1.0 - (2 * height / (height * height + a * a)));
        g_inf = sqrt(Sigma * Delta / rho_squared);
      }
    }
  }
  // printf("energy shift ring source - observer (ONLY e^-nu NOW): %e, h, r, a = %e, %e, %e\n", g_inf, height, r_ring, a);
  return g_inf;
}

double calc_cos_theta_primary_source(double height, double r_ring, double a) {
  double cosine_theta;
  if (r_ring == 0.0) {
    cosine_theta = 1.0; // on the rotational axis
  }
  else if (a == 0.0){ // special case of zero spin
    cosine_theta = height / sqrt(height * height + r_ring * r_ring);
  } else { // general case, the solution of equation r_ring = sqrt(h^2 * tan(theta)^2 + a^2 * sin(theta)^2) for cos(theta)
    cosine_theta = sqrt((sqrt(4.0 * a * a * height * height + (height * height + r_ring * r_ring - a * a) *
        (height * height + r_ring * r_ring - a * a)) + a * a - height * height - r_ring * r_ring) / (2.0 * a * a));
  }
  return cosine_theta;
}

double calc_spherical_radius_primary_source(double height, double r_ring, double a) {
  double cosine_theta = calc_cos_theta_primary_source(height, r_ring, a);
  double r = height / cosine_theta; // the spherical radius by definition
  return r;
}

// double convert_sph_radius_costheta_to_h_rring()

double calc_extent_from_sph_radius_costheta(double sph_radius, double costheta, double a) {
  return sqrt(1.0 - costheta * costheta) * sqrt(sph_radius * sph_radius + a * a);
}

double calc_extent_from_sph_radius_primary_source(double height, double sph_radius, double a) {
  if (abs(height - sph_radius) < 1e-4){
    return 0.0; // it is lamp post
  }
  assert(height <= sph_radius && ("Height cannot be more than spherical radius"));
  return sqrt((sph_radius * sph_radius - height * height) * (1.0 + a * a / (sph_radius * sph_radius)));
}


/* get RMS (ISCO) for the Kerr Case */
double kerr_rms(double a) {
  //	 accounts for negative spin
  double sign = 1.0;
  if (a < 0) {
    sign = -1.0;
  }

  double Z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0));
  double Z2 = sqrt((3.0 * a * a) + (Z1 * Z1));

  return 3.0 + Z2 - sign * sqrt((3.0 - Z1) * (3.0 + Z1 + (2 * Z2)));
}

/* get the rplus value (size of the black hole event horizon */
double kerr_rplus(double a) {
  return 1.0 + sqrt(1.0 - a * a);
}

/** calculate the doppler factor for a moving primary source **/
double doppler_factor(double del, double bet) {
  return sqrt(1.0 - bet * bet) / (1.0 + bet * cos(del));
}

/** similar to gi_potential_lp below, but for ring source **/
double gi_potential_ext(double r, double a, double h, double x) {
  // as a temporary solution, exploit previously written function applicable to this case
  // in general I need to take it from Rellp_Extended, as there we have proper g_sd source-disk
  double g_inf_s = calc_g_inf_ext(h, x, a); // = E_inf / E_s
  double g_inf_d = calc_g_inf_ext(0.0, r, a); // = E_inf / E_d
  return g_inf_s / g_inf_d * gamma_correction(r, a); // = E_d / E_s
}

/** calculates g = E/E_i in the lamp post geometry (see, e.g., 27 in Dauser et al., 2013, MNRAS) **/
double gi_potential_lp(double r, double a, double h, double bet, double del) {

  /** ! calculates g = E/E_i in the lamp post geometry
    ! (see, e.g., page 48, Diploma Thesis, Thomas Dauser) **/
  double ut_d = ((r * sqrt(r) + a) / (sqrt(r) * sqrt(r * r - 3 * r + 2 * a * sqrt(r))));
  double ut_h = sqrt((h * h + a * a) / (h * h - 2 * h + a * a));

  double gi = ut_d / ut_h;

  // check if we need to calculate the additional factor for the velocity
  if (fabs(bet) < 1e-6) {
    return gi;
  }

  double gam = 1.0 / sqrt(1.0 - bet * bet);

  // get the sign for the equation
  double sign = 1.0;
  if (del > M_PI / 2) {
    sign = -1.0;
  }

  // ** Adam, 28.7.2021: **
  // So given the expression for Carter’s q:
  //q^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h - a^2,
  //it is clear that
  //q^2 + a^2 = \sin^2\delta (h^2+a^2)^2/\Delta_h
  //Therefore:
  //(h^2+a^2)^2 - \Delta_h (q^2+a^2) = (h^2+a^2)^2 - \sin^2\delta (h^2+a^2)^2
  //= (h^2+a^2)^2 \cos^2\delta
  //
  //Therefore:
  //\sqrt{ (h^2+a^2)^2 - \Delta_h (q^2+a^2) } / (h^2+a^2) = \cos\delta
  //
  //Therefore your equation (27) becomes:
  //glp = glp(beta=0) / { \gamma [ 1 -/+ \beta\cos\delta ] }

  double delta_eq = h * h - 2 * h + a * a;
  double q2 = (pow(sin(del), 2)) * (pow((h * h + a * a), 2) / delta_eq) - a * a;

  double beta_fac = sqrt(pow((h * h + a * a), 2) - delta_eq * (q2 + a * a));
  beta_fac = gam * (1.0 + sign * beta_fac / (h * h + a * a) * bet);

  return gi / beta_fac;
}

/** similar to gi_potential_lp below, but for ring source **/
double gi_potential_ext(const relParam *param, double r) {
  // temporary older solution
  // double g_inf_s = calc_g_inf_ext(param->height, param->r_ring, param->a); // = E_inf / E_s
  // double g_inf_d = calc_g_inf_ext(0.0, r, param->a); // = E_inf / E_d
  //return g_inf_s / g_inf_d * gamma_correction(r, param->a); // = E_d / E_s
  double g_sd = get_energyboost_ext_source_disk(param, r, 1.0); // it is just energy shift here, so no Gamma
  // printf("Do we agree at LP? g old: %.6f, new %.6f \n", g_inf_s / g_inf_d * gamma_correction(r, param->a), g_sd);
  return g_sd;
}

/**
 * @brief calculate the doppler factor from source to the observer
 * we currently assume that delta_obs and incl are the same (so no GR light-bending
 * is accounted for).
 * @param incl [in rad]
 * @param beta [in v/c]
 * @return
 */
double doppler_factor_source_obs(const relParam *rel_param) { //update?
  // note that the angles incl and delta_obs are defined the opposite way
  double delta_obs = M_PI - rel_param->incl;
  assert(delta_obs > M_PI / 2);
  assert(delta_obs < M_PI);

  // glp = D*glp(beta=0) = glp(beta=0) / { \gamma [ 1 -/+ \beta\cos\delta ] }
  // Adam: we use the minus sign for δ > π/2 and the plus sign for δ < π/2 to get
  // -> as always δ > π/2 in this case, we directly can use the minus sign
  double doppler = doppler_factor(delta_obs, rel_param->beta);
  // assert(doppler >= 1-1e-8);

  return doppler;
}

/**
 * @brief calculate the energy shift from the LP to the observer, also taking
 * a potential velocity of the source into account
 *
 * note that this routine uses  the function "doppler_factor_source_obs", which
 * sets δ_obs = incl, which is not fully true in GR, however, within 1% accuracy
 */
double energy_shift_source_obs(const relParam *rel_param) {

  double g_inf_0 = 1;
  // if we do NOT have the LP geometry, the energy shift is just 1
  // (as without geometry we can not define an energy shift)
  if (rel_param == nullptr || rel_param->emis_type != EMIS_TYPE_LP) {
    return g_inf_0;
  }

  if (rel_param->model_type == MOD_TYPE_RELXILLLPEXT) {
    // g_inf_0 = calc_g_inf_ext(rel_param->height, rel_param->r_ring, rel_param->a);
    g_inf_0 = get_energyboost_ext_source_obs(rel_param, 1.0);
  } else {
    g_inf_0 = calc_g_inf(rel_param->height, rel_param->a);
  }

  if (rel_param->beta < 1e-4) {
    return g_inf_0;
  } else {
    return g_inf_0 * doppler_factor_source_obs(rel_param);
  }

}

/**
 * @brief calculate the energy shift from the LP to the disk zone with mean radius radius_disk of photon trajectories
 * emitted under del_emit (in the frame of the non-moving source); del_emit is only important for beta>0
 * @details AN: currently used in IonGradient.cpp
 * The connection for which del_emit which radius radius_disk is hit, is stored in the LP table structure
 */
double energy_shift_source_disk(const relParam *rel_param, double radius_disk, double del_emit) {

  // if we do NOT have the LP geometry, the energy shift is just 1
  // (as without geometry we can not define an energy shift)
  if (rel_param == nullptr || rel_param->emis_type != EMIS_TYPE_LP) {
    return 1;
  }
  if (rel_param->model_type == MOD_TYPE_RELXILLLPEXT) {
    double gshift = get_energyboost_ext_source_disk(rel_param, radius_disk, 1.0);
    // return gi_potential_ext(radius_disk, rel_param->a, rel_param->height, rel_param->r_ring);
    return gshift;
  }
  else{
    return gi_potential_lp(radius_disk, rel_param->a, rel_param->height, rel_param->beta, del_emit);
  }
}

/**
 * @brief calculate the geometric/Newtonian emissivity for the LP geometry
 * @details the emissivity profiles are normalized such that at large height and
 * radii the GR profiles converge towards this defintion
 * @param h
 * @param r
 * @return
 */
double calc_lp_emissivity_newton(double h, double r) {

  double emis = pow(1. / (pow((r / h), 2) + 1), (3. / 2));
  emis /= (2 * M_PI * h * h);

  return emis;
}

/**
 * @brief calculate the fluxboost from the (moving) primary source to the accretion disk
 * at radius rad, following g^Gamma*Doppler_factor^2 (see Dauser et al., 2013)
 * @details AN: used in Rellp to normalize flux of the LP model only. For ext models we do this during interpolation
 * Which is a bit inconsistent, so TODO (maybe) unify boost calculation for lp and ext models
 * @param rad
 * @param del_emit
 * @param a
 * @param height
 * @param gamma
 * @param beta
 * @return
 */
double calc_fluxboost_source_disk(double rad, double del_emit, double a, double height, double gamma, double beta) {

  // then energy shift factor gi^gamma (see Dauser et al., 2013)
  double boost = pow(gi_potential_lp(rad, a, height, beta, del_emit), gamma);

  // take the beaming of the jet into account (see Dauser et al., 2013)
  if (beta > 1e-6) {
    boost *= pow(doppler_factor(del_emit, beta), 2);
  }
  return boost;
}

double kerr_to_flat_area_ratio(double r, double a) {
  // the value below is obtained directly from ratio A_Kerr / A_flat by all possible simplifications
  return sqrt( (r * r * r + a * a * r + 2.0 * a * a) /
                  (r * r * r + a * a * r - 2.0 * r * r) );
}

double gamma_correction(double r, double a){
  // the value is obtained from the "Bardeen"-like form of the Lorentz factor by a few simplifications
  return (r * sqrt(r) + abs(a)) * sqrt( (r * r - 2.0 * r + a * a) / (r * r +
              2.0 * a * sqrt(r) - 3.0 * r) / (r * r * r + a * a * r + 2.0 * a * a) );
}
