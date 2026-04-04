#include "CHICEARTH.h"
#include <cmath>
#include <cassert>
 
// =================================================== //
// ====== Class for Earth neutrino propagation ======= //
// =================================================== //
 
CHICEARTH::CHICEARTH(std::string_view mode,
                     double theta_12, double theta_23, double theta_13,
                     double delta_cp, double dm2_21, double dm2_31,
                     std::string_view model)

    : chic(std::make_unique<CHIC>(mode, theta_12, theta_23, theta_13,
                                  delta_cp, dm2_21, dm2_31))
    , Earth(&get_earth_model(model))
{
  assert(Earth->Nlayers > 0 && "Earth model has no layers");
 
  // Pre-allocate layer cache and track buffers — avoids repeated heap use later
  Layers.resize(Earth->Nlayers);
  tracks.resize(Earth->Nlayers, 0.0);
  radii2.resize(Earth->Nlayers);
  for (int i = 0; i < Earth->Nlayers; ++i)
    radii2[i] = Earth->radii[i] * Earth->radii[i];
 
  // Initialise working matrices to identity / zero
  J.setIdentity();
  J_layer.setZero();
}

// Computes the eigenvalues and Hamiltonian for each layer at a given energy
void CHICEARTH::_compute_hamiltonians() {
  if (need_update) {
    for (int i = 0; i < Earth->Nlayers; i++) {
      // Update CHIC with layer parameters
      chic->update_density(Earth->density[i]);
      chic->update_Ye(Earth->Ye[i]);
      chic->compute_hamiltonians(E0);
      // Store relevant energy-dependent quantities for the layer
      Layers[i].lambdas = chic->get_eigenvalues();
      Layers[i].prod_lambdas = chic->get_prod_eigenvalues();
      Layers[i].diff_lambdas = chic->get_diff_eigenvalues();
      Layers[i].Hs = chic->get_hamiltonian();
      Layers[i].Hs2 = chic->get_hamiltonian_squared();
    }
    need_update = false;
  }
}

void CHICEARTH::_build_track() {
  deepest = Earth->Nlayers - 1; // starts in outermost shell
  const double Rsin2 = R_EARTH * R_EARTH * std::max(1.0 - cos_zenith0 * cos_zenith0, 0.0);
  
  for (int i = 0; i < Earth->Nlayers; ++i) {  // outer → inner
    const double seg2 = radii2[i] - Rsin2;
    if (seg2 > 0.0) {
      tracks[i] = std::sqrt(seg2);
      if (deepest > i) deepest = i;  // updated to the innermost hit layer
    } else {
      tracks[i] = 0.0;
    }
  }
}


void CHICEARTH::_layer_amplitude(double L, const Layer& lay) {
  // Exponentials
  iL = std::complex<double>(0.0, -L * OptConstants::BASELINE_FACTOR);
  exp_eigenvals = lay.diff_lambdas.array() * (iL * lay.lambdas).array().exp();
 
  // Cayley-Hamilton coefficients
  const std::complex<double> c2 = exp_eigenvals.sum();
  const std::complex<double> c1 = lay.lambdas.cast<std::complex<double>>().dot(exp_eigenvals);
  // c0 is the dot with prod_lambdas (diagonal / identity coefficient)
  const std::complex<double> c0 = lay.prod_lambdas.cast<std::complex<double>>().dot(exp_eigenvals);
 
  J_layer.noalias() = c2 * lay.Hs2;
  J_layer.noalias() += c1 * lay.Hs;
  J_layer.diagonal().array() += c0;

  // Sandwich product for combined amplitude
  J.applyOnTheLeft(J_layer);
  J.applyOnTheRight(J_layer);
}

void CHICEARTH::_amplitude() {
  J.setIdentity();
  // start with the deepest layer
  _layer_amplitude(2 * tracks[deepest], Layers[deepest]);
  // loop over the rest of shallower layers
  for (int i = deepest + 1; i < Earth->Nlayers; i++) {
    _layer_amplitude(std::abs(tracks[i] - tracks[i - 1]), Layers[i]);
  }
}

Eigen::Matrix3d CHICEARTH::compute_oscillations(double E, double cos_zenith) {
  if (cos_zenith >= 0.0)
    return Eigen::Matrix3d::Identity(); // to be removed when h is applied

  // Energy dependence
  if (E != E0) {
    E0 = E;
    need_update = true;
  }
  _compute_hamiltonians();

  // Build track
  if (cos_zenith != cos_zenith0) {
    cos_zenith0 = cos_zenith;
    _build_track();
  }

  // Baseline dependence and amplitude
  _amplitude();

  return J.cwiseAbs2();
}
