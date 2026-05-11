#include "CHIC_EARTH.h"
#include "earth.h"
#include <cmath>
#include <cassert>

#include <iostream>
 
// =================================================== //
// ====== Class for Earth neutrino propagation ======= //
// =================================================== //
 
CHICEARTH::CHICEARTH(std::string_view mode,
                     double theta_12, double theta_23, double theta_13,
                     double delta_cp, double dm2_21, double dm2_31,
                     std::string_view earth_model, double detector_depth)

    : chic(std::make_unique<CHIC>(mode, theta_12, theta_23, theta_13,
                                  delta_cp, dm2_21, dm2_31))
    , Earth(&get_earth_model(earth_model))
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
}

void CHICEARTH::_build_track() {
  deepest = Earth->Nlayers - 1; // starts in outermost shell
  const double Rsin2 = R2_EARTH * std::max(1.0 - cos_zenith0 * cos_zenith0, 0.0);
  
  for (int i = 0; i < Earth->Nlayers; ++i) {  // outer → inner
    const double seg2 = radii2[i] - Rsin2;
    if (seg2 > 0.0) {
      tracks[i] = std::sqrt(seg2);
      if (i < deepest) deepest = i;
    }
    else tracks[i] = 0.0;
  }
}

void CHICEARTH::_layer_amplitude(double L, const Layer& lay, const bool deepest) {
  // Exponentials
  iL = std::complex<double>(0.0, -L * OptConstants::BASELINE_FACTOR);
  exp_eigenvals = lay.diff_lambdas.array() * (iL * lay.lambdas).array().exp();
  
  J_layer.noalias() = exp_eigenvals.sum() * lay.Hs2;
  J_layer.noalias() += lay.lambdas.dot(exp_eigenvals) * lay.Hs;
  J_layer.diagonal().array() += lay.prod_lambdas.dot(exp_eigenvals);

  // Sandwich product for combined amplitude
  if (deepest) J = J_layer;
  else J = J_layer * J * J_layer;
  // J.applyOnTheLeft(J_layer);
  // J.applyOnTheRight(J_layer);
}

void CHICEARTH::_amplitude() {
  // amplitude at deepest layer
  _layer_amplitude(2 * tracks[deepest], Layers[deepest], true);
  // loop over the rest of shallower layers
  for (int i = deepest + 1; i < Earth->Nlayers; i++) {
    _layer_amplitude(std::abs(tracks[i] - tracks[i - 1]), Layers[i]);
  }
}

Eigen::Matrix3d CHICEARTH::compute_oscillations(double E, double cos_zenith, double h) {
  if (cos_zenith >= 0.0)
    return Eigen::Matrix3d::Identity(); // to be removed when h is applied

  // Build track
  if (cos_zenith != cos_zenith0) {
    cos_zenith0 = cos_zenith;
    _build_track();
  }

  // Energy dependence
  if (E != E0 || update_vacuum) {
    E0 = E;
    _compute_hamiltonians();
  }

  // Baseline dependence and amplitude
  _amplitude();
 return J.cwiseAbs2();
}


// ====================================================================== \\
// ==== Class for probabilities and derivatives of Earth propagation ==== \\
// ====================================================================== \\
 
CHICDIFFEARTH::CHICDIFFEARTH(std::string_view mode,
                     double theta_12, double theta_23, double theta_13,
                     double delta_cp, double dm2_21, double dm2_31,
                     std::string_view earth_model, double detector_depth)

    : dchic(std::make_unique<CHICDIFF>(mode, theta_12, theta_23, theta_13,
                                  delta_cp, dm2_21, dm2_31))
    , Earth(&get_earth_model(earth_model))
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
void CHICDIFFEARTH::_compute_hamiltonians_and_anticommutators() {
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
}

void CHICDIFFEARTH::_build_track() {
  deepest = Earth->Nlayers - 1; // starts in outermost shell
  const double Rsin2 = R2_EARTH * std::max(1.0 - cos_zenith0 * cos_zenith0, 0.0);
  
  for (int i = 0; i < Earth->Nlayers; ++i) {  // outer → inner
    const double seg2 = radii2[i] - Rsin2;
    if (seg2 > 0.0) {
      tracks[i] = std::sqrt(seg2);
      if (i < deepest) deepest = i;
    }
    else tracks[i] = 0.0;
  }
}

void CHICDIFFEARTH::_layer_amplitude_diff(double L, const Layer& lay, const bool deepest) {
  // Integral matrix (symmetric)
    I_ijk.diagonal() = iL * lay.diff_lambdas.cwiseProduct(exp_eigenvals);
    I_ijk(0, 1) = I_ijk(1, 0) = (lay.diff_lambdas[0]*exp_eigenvals[1] - lay.diff_lambdas[1]*exp_eigenvals[0])
                / (lay.lambdas[1] - lay.lambdas[0]);
    I_ijk(0, 2) = I_ijk(2, 0) = (lay.diff_lambdas[2]*exp_eigenvals[0] - lay.diff_lambdas[0]*exp_eigenvals[2])
                / (lay.lambdas[0] - lay.lambdas[2]);
    I_ijk(1, 2) = I_ijk(2, 1) = (lay.diff_lambdas[1]*exp_eigenvals[2] - lay.diff_lambdas[2]*exp_eigenvals[1])
                / (lay.lambdas[2] - lay.lambdas[1]);

    Eigen::Vector3cd I_prod = I_ijk * prod_lambdas;

    cdJ_diag[0] = prod_lambdas.dot(I_prod);    // dH

    I_prod = I_ijk * lambdas;
    cdJ_diag[1] = lambdas.dot(I_prod);       // H dH H
    cdJ_off[0] = prod_lambdas.dot(I_prod);   // dH H

    I_prod = I_ijk * unit;
    cdJ_diag[2] = unit.dot(I_prod);            // H2 dH H2
    cdJ_off[1] = prod_lambdas.dot(I_prod);     // dH H2
    cdJ_off[2] = lambdas.dot(I_prod);          // H dH H2

    dJ_layer = cdJ_diag[0] * dlay.dHs;
    dJ_layer.noalias()+= cdJ_diag[1] * dlay.Hs_dHs_Hs;
    dJ_layer.noalias()+= cdJ_diag[2] * dlay.Hs2_dHs_Hs2;
    dJ_layer.noalias()+= cdJ_off[0] * dlay.comm_dHH;
    dJ_layer.noalias()+= cdJ_off[1] * dlay.comm_dHH2;
    dJ_layer.noalias()+= cdJ_off[2] * dlay.Hs_comm_dHH_Hs;

  // Sandwich product for combined amplitude and the corresponding for the derivative
  if (deepest) {
    dJ = dJ_layer
    J = J_layer;
  }
  else {
    dJ = J_layer * dJ * J_layer + J_layer * J * dJ_layer + dJ_layer * J * dJ_layer;
    J = J_layer * J * J_layer;
  }
  // J.applyOnTheLeft(J_layer);
  // J.applyOnTheRight(J_layer);
}


void CHICDIFFEARTH::_amplitude_and_diff() {
  _layer_amplitude_and_diff(2 * tracks[deepest], Layers[deepest], true);
  // loop over the rest of shallower layers
  for (int i = deepest + 1; i < Earth->Nlayers; i++) {
    std::cout<<"Density: "<<Earth->density[i]<<std::endl;
    _layer_amplitude(std::abs(tracks[i] - tracks[i - 1]), Layers[i]);
    _layer_amplitude_and_diff(std::abs(tracks[i] - tracks[i - 1]), Layers[i]);
  }
}

Eigen::Matrix3d CHICDIFFEARTH::compute_oscillations(double E, double cos_zenith, double h) {
  if (cos_zenith >= 0.0)
    return Eigen::Matrix3d::Identity(); // to be removed when h is applied

  // Build track
  if (cos_zenith != cos_zenith0) {
    cos_zenith0 = cos_zenith;
    _build_track();
  }

  // Energy dependence
  if (E != E0 || update_vacuum) {
    E0 = E;
    _compute_hamiltonians_and_anticommutators();
  }

  // Baseline dependence and amplitude
  _amplitude();
 return J.cwiseAbs2();
}

Eigen::Matrix3d CHICDIFFEARTH::compute_oscillations_derivatives(std::string_view param, double E, double cos_zenith, double h) {
  if (cos_zenith >= 0.0)
    return Eigen::Matrix3d::Identity(); // to be removed when h is applied

  // Build track
  if (cos_zenith != cos_zenith0) {
    cos_zenith0 = cos_zenith;
    _build_track();
  }

  // Energy dependence
  if (E != E0 || update_vacuum) {
    E0 = E;
    _compute_hamiltonians_and_anticommutators();
  }

  // Baseline dependence and amplitude
  _amplitude_and_diff();

  return 2.0 * dJ.cwiseProduct(J.conjugate()).real();
}



/*
To compute derivative, we need to compute the amplitude.
Make up your mind first if you want to get the derivatives.
So we can implement the calculation of J and dJ at the same time and save time.
Maybe split in files? Not really. 
Actually we don't need full dJ, but 
*/
