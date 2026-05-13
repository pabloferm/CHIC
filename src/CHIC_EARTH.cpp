#include "CHIC_EARTH.h"
#include "CHIC.h"
#include "earth.h"
#include <cassert>
 
// =================================================== //
// ====== Class for Earth neutrino propagation ======= //
// =================================================== //
 
CHICEARTH::CHICEARTH(std::unique_ptr<CHIC> chic_ptr,
                     std::string_view earth_model, double detector_depth)
    : chic(std::move(chic_ptr))
    , Earth(&get_earth_model(earth_model))
{
  assert(Earth->Nlayers > 0 && "Earth model has no layers");

  Layers.resize(Earth->Nlayers);
  tracks.resize(Earth->Nlayers, 0.0);
  radii2.resize(Earth->Nlayers);
  for (int i = 0; i < Earth->Nlayers; ++i)
    radii2[i] = Earth->radii[i] * Earth->radii[i];

  J.setIdentity();
  J_layer.setZero();
}

// --- Public constructor to be compatible with CHIC and CHCDIFF
CHICEARTH::CHICEARTH(std::string_view mode,
                     double theta_12, double theta_23, double theta_13,
                     double delta_cp, double dm2_21, double dm2_31,
                     std::string_view earth_model, double detector_depth)
    : CHICEARTH(std::make_unique<CHIC>(mode, theta_12, theta_23, theta_13,
                                       delta_cp, dm2_21, dm2_31),
                earth_model, detector_depth)
{}

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
  update_vacuum = false;
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

void CHICEARTH::_layer_amplitude(int i_layer, const bool deepest) {
  // Exponentials
  
  if (deepest) iL = std::complex<double>(0.0, -2*tracks[i_layer] * OptConstants::BASELINE_FACTOR);
  else iL = std::complex<double>(0.0, -std::abs(tracks[i_layer] - tracks[i_layer - 1]) * OptConstants::BASELINE_FACTOR);

  exp_eigenvals = Layers[i_layer].diff_lambdas.array() * (iL * Layers[i_layer].lambdas).array().exp();
  
  J_layer.noalias() = exp_eigenvals.sum() * Layers[i_layer].Hs2;
  J_layer.noalias() += Layers[i_layer].lambdas.dot(exp_eigenvals) * Layers[i_layer].Hs;
  J_layer.diagonal().array() += Layers[i_layer].prod_lambdas.dot(exp_eigenvals);
}


void CHICEARTH::_update_amplitude(bool deepest) {
  // Sandwich product for combined amplitude
  if (deepest) J = J_layer;
  else J = J_layer * J * J_layer;
}

void CHICEARTH::_amplitude() {
  // amplitude at deepest layer
  _layer_amplitude(deepest, true);
  _update_amplitude(true);
  // loop over the rest of shallower layers
  for (int i = deepest + 1; i < Earth->Nlayers; i++) {
    _layer_amplitude(i);
    _update_amplitude();
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


// ======================================================================
// ==== Class for probabilities and derivatives of Earth propagation ====
// ======================================================================
 
CHICEARTHDIFF::CHICEARTHDIFF(std::string_view mode,
                             double theta_12, double theta_23, double theta_13,
                             double delta_cp, double dm2_21, double dm2_31,
                             std::string_view earth_model, double detector_depth)
    : CHICEARTH(std::make_unique<CHICDIFF>(mode, theta_12, theta_23, theta_13,
                                           delta_cp, dm2_21, dm2_31),
                earth_model, detector_depth)
{
  dLayers.resize(Earth->Nlayers);
  dJ.setZero();
  dJ_layer.setZero();
}

// Computes eigenvalues, Hamiltonians, and anticommutators for each layer
void CHICEARTHDIFF::_compute_hamiltonians_and_anticommutators() {
  auto* chic_diff = static_cast<CHICDIFF*>(chic.get());
  for (int i = 0; i < Earth->Nlayers; i++) {
    chic_diff->update_density(Earth->density[i]);
    chic_diff->update_Ye(Earth->Ye[i]);
    chic_diff->compute_hamiltonians_and_derivatives(param0, E0);
    Layers[i].lambdas           = chic_diff->get_eigenvalues();
    Layers[i].prod_lambdas      = chic_diff->get_prod_eigenvalues();
    Layers[i].diff_lambdas      = chic_diff->get_diff_eigenvalues();
    Layers[i].Hs                = chic_diff->get_hamiltonian();
    Layers[i].Hs2               = chic_diff->get_hamiltonian_squared();
    dLayers[i].dHs              = chic_diff->get_diff_hamiltonian();
    dLayers[i].Hs_dHs_Hs        = Layers[i].Hs * dLayers[i].dHs * Layers[i].Hs;
    dLayers[i].Hs2_dHs_Hs2      = Layers[i].Hs2 * dLayers[i].dHs * Layers[i].Hs2;
    dLayers[i].comm_dHsHs       = chic_diff->get_comm_dHsHs();
    dLayers[i].comm_dHsHs2      = chic_diff->get_comm_dHsHs2();
    dLayers[i].Hs_comm_dHsHs_Hs = Layers[i].Hs * dLayers[i].comm_dHsHs * Layers[i].Hs;
  }
  update_vacuum = false;
  update_param = false;
}


void CHICEARTHDIFF::_layer_amplitude_diff(int i_layer, const bool deepest) {
  // Integral matrix (symmetric)
    I_ijk.diagonal() = iL * Layers[i_layer].diff_lambdas.cwiseProduct(exp_eigenvals);
    I_ijk(0, 1) = I_ijk(1, 0) = (Layers[i_layer].diff_lambdas[0]*exp_eigenvals[1] - Layers[i_layer].diff_lambdas[1]*exp_eigenvals[0])
                / (Layers[i_layer].lambdas[1] - Layers[i_layer].lambdas[0]);
    I_ijk(0, 2) = I_ijk(2, 0) = (Layers[i_layer].diff_lambdas[2]*exp_eigenvals[0] - Layers[i_layer].diff_lambdas[0]*exp_eigenvals[2])
                / (Layers[i_layer].lambdas[0] - Layers[i_layer].lambdas[2]);
    I_ijk(1, 2) = I_ijk(2, 1) = (Layers[i_layer].diff_lambdas[1]*exp_eigenvals[2] - Layers[i_layer].diff_lambdas[2]*exp_eigenvals[1])
                / (Layers[i_layer].lambdas[2] - Layers[i_layer].lambdas[1]);

    Eigen::Vector3cd I_prod = I_ijk * Layers[i_layer].prod_lambdas;

    cdJ_diag[0] = Layers[i_layer].prod_lambdas.dot(I_prod);    // dH

    I_prod = I_ijk * Layers[i_layer].lambdas;
    cdJ_diag[1] = Layers[i_layer].lambdas.dot(I_prod);       // H dH H
    cdJ_off[0] = Layers[i_layer].prod_lambdas.dot(I_prod);   // dH H

    I_prod = I_ijk * Eigen::Vector3cd::Ones();
    cdJ_diag[2] = Eigen::Vector3cd::Ones().dot(I_prod);            // H2 dH H2
    cdJ_off[1] = Layers[i_layer].prod_lambdas.dot(I_prod);     // dH H2
    cdJ_off[2] = Layers[i_layer].lambdas.dot(I_prod);          // H dH H2

    dJ_layer = cdJ_diag[0] * dLayers[i_layer].dHs;
    dJ_layer.noalias()+= cdJ_diag[1] * dLayers[i_layer].Hs_dHs_Hs;
    dJ_layer.noalias()+= cdJ_diag[2] * dLayers[i_layer].Hs2_dHs_Hs2;
    dJ_layer.noalias()+= cdJ_off[0] * dLayers[i_layer].comm_dHsHs;
    dJ_layer.noalias()+= cdJ_off[1] * dLayers[i_layer].comm_dHsHs2;
    dJ_layer.noalias()+= cdJ_off[2] * dLayers[i_layer].Hs_comm_dHsHs_Hs;

  // Sandwich product for combined amplitude and the corresponding for the derivative
  if (deepest)  dJ = dJ_layer;
  else dJ = J_layer * dJ * J_layer + J_layer * J * dJ_layer + dJ_layer * J * J_layer;
}

void CHICEARTHDIFF::_amplitude_and_diff() {
  _layer_amplitude(deepest, true);
  _layer_amplitude_diff(deepest, true);
  _update_amplitude(true);
  // loop over the rest of shallower layers
  for (int i = deepest + 1; i < Earth->Nlayers; i++) {
    _layer_amplitude(i);
    _layer_amplitude_diff(i);
    _update_amplitude();
  }
}

Eigen::Matrix3d CHICEARTHDIFF::compute_oscillations(double E, double cos_zenith, double h) {
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
    update_param = true;
  }

  // Baseline dependence and amplitude
  _amplitude();
 return J.cwiseAbs2();
}

Eigen::Matrix3d CHICEARTHDIFF::compute_oscillations_derivatives(std::string_view param, double E, double cos_zenith, double h) {
  if (param != param0) {
    param0        = std::string(param);
    update_param = true;
  }

  if (cos_zenith >= 0.0)
    return Eigen::Matrix3d::Identity(); // to be removed when h is applied

  // Build track
  if (cos_zenith != cos_zenith0) {
    cos_zenith0 = cos_zenith;
    _build_track();
  }

  // Energy dependence
  if (E != E0 || update_vacuum || update_param) {
    E0 = E;
    _compute_hamiltonians_and_anticommutators();
  }

  // Baseline dependence and amplitude
  _amplitude_and_diff();

  return 2.0 * dJ.cwiseProduct(J.conjugate()).real();
}
