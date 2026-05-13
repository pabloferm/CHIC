#ifndef CHICEARTH_H
#define CHICEARTH_H
 
#include "CHIC.h"
#include "earth.h"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

// ===================================================
// ====== Class for Earth neutrino propagation =======
// ===================================================

class CHICEARTH {

public:
  explicit CHICEARTH(
      std::string_view mode    = "neutrino",
      double theta_12          = 0.5836381018669037,
      double theta_23          = 0.8587019919812102,
      double theta_13          = 0.14957471689591406,
      double delta_cp          = 4.084070449666731,
      double dm2_21            = 7.42e-5,
      double dm2_31            = 2.51e-3,
      std::string_view model   = "PREM10",
      double detector_depth    = 0.0);

  // ---- Parameter update setters ----
  inline void update_dcp  (double delta_cp)  { chic->update_dcp(delta_cp);   update_vacuum = true; }
  inline void update_th23 (double theta_23)  { chic->update_th23(theta_23);  update_vacuum = true; }
  inline void update_th12 (double theta_12)  { chic->update_th12(theta_12);  update_vacuum = true; }
  inline void update_th13 (double theta_13)  { chic->update_th13(theta_13);  update_vacuum = true; }
  inline void update_dm221(double dm_2_21)   { chic->update_dm221(dm_2_21);  update_vacuum = true; }
  inline void update_dm231(double dm_2_31)   { chic->update_dm231(dm_2_31);  update_vacuum = true; }

  virtual Eigen::Matrix3d compute_oscillations(double E, double cos_zenith, double h = 0.0);

  [[nodiscard]] const Eigen::Matrix3cd& get_amplitude() const noexcept { return J; }

  virtual ~CHICEARTH() = default;

protected:

  // Injection constructor — used by derived classes to supply a CHIC subclass
  CHICEARTH(std::unique_ptr<CHIC> chic_ptr,
            std::string_view earth_model, double detector_depth);

  // ---- Per-layer cached quantities (energy-dependent) ----
  struct Layer {
    Eigen::Matrix3cd Hs;            // traceless H
    Eigen::Matrix3cd Hs2;           // traceless H²
    Eigen::Vector3d  lambdas;       // eigenvalues
    Eigen::Vector3d  diff_lambdas;  // eigenvalue differences
    Eigen::Vector3d  prod_lambdas;  // pairwise products of eigenvalues
  };

  void _compute_hamiltonians();
  void _build_track();
  void _amplitude();
  // Two-phase layer update — split so derived class can insert dJ update between them:
  void _compute_J_layer(int i_layer, bool is_deepest = false); // sets iL, exp_eigenvals, J_layer
  void _update_J       (bool is_deepest = false);              // applies J = J_layer * J * J_layer
  void _layer_amplitude(int i_layer, bool is_deepest = false); // calls both (convenience for base)

  // ---- Data members ----
  std::unique_ptr<CHIC> chic;
  const EarthModel*     Earth = nullptr;

  std::vector<Layer>  Layers;
  std::vector<double> tracks, radii2;

  Eigen::Matrix3cd J;
  Eigen::Matrix3cd J_layer;
  Eigen::Vector3cd exp_eigenvals;

  bool   update_vacuum = true;
  int    deepest       = 0;
  double cos_zenith0{-2.0};  // out of valid range [-1, 0)
  double E0{-1.0};           // out of valid range
  const double R2_EARTH = R_EARTH * R_EARTH;
  std::complex<double> iL;
};


// ===================================================
// == Derived class: probabilities + derivatives =====
// ===================================================

class CHICEARTHDIFF : public CHICEARTH {
public:
  explicit CHICEARTHDIFF(
      std::string_view mode    = "neutrino",
      double theta_12          = 0.5836381018669037,
      double theta_23          = 0.8587019919812102,
      double theta_13          = 0.14957471689591406,
      double delta_cp          = 4.084070449666731,
      double dm2_21            = 7.42e-5,
      double dm2_31            = 2.51e-3,
      std::string_view model   = "PREM10",
      double detector_depth    = 0.0);

  Eigen::Matrix3d compute_oscillations(double E, double cos_zenith, double h = 0.0) override;
  Eigen::Matrix3d compute_oscillations_derivatives(std::string_view param, double E, double cos_zenith, double h = 0.0);

private:

  struct dLayer {
    Eigen::Matrix3cd dHs;
    Eigen::Matrix3cd Hs_dHs_Hs;
    Eigen::Matrix3cd Hs2_dHs_Hs2;
    Eigen::Matrix3cd comm_dHsHs;         // {dH, H}  = dH*H + H*dH
    Eigen::Matrix3cd comm_dHsHs2;        // {dH, H²} = dH*H² + H²*dH
    Eigen::Matrix3cd Hs_comm_dHsHs_Hs;  // H * {dH,H} * H
  };

  void _compute_hamiltonians_and_anticommutators();
  void _amplitude_and_diff();
  // Order in _amplitude_and_diff: _compute_J_layer → _layer_amplitude_diff → _update_J.
  // dJ needs J_old (before sandwich) and J_layer, so it must come before _update_J.
  void _layer_amplitude_diff(int i_layer, bool is_deepest = false);

  std::vector<dLayer> dLayers;

  // FIX: track current param so we detect changes and recompute dLayers
  std::string param0{"none"};
  bool        param_changed = true;

  Eigen::Matrix3cd dJ;
  Eigen::Matrix3cd dJ_layer;
  Eigen::Matrix3cd I_ijk;
  std::complex<double> cdJ_diag[3], cdJ_off[3];
};

#endif // CHICEARTH_H
