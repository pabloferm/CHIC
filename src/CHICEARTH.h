#ifndef CHICEARTH_H
#define CHICEARTH_H
 
#include "CHIC.h"
#include "earth.h"
#include <Eigen/Dense>
#include <memory>   // std::unique_ptr
#include <vector>

// =================================================== \\
// ====== Class for Earth neutrino propagation ======= \\
// =================================================== \\

class CHICEARTH {

public:
  // Constructor
  explicit CHICEARTH(
      std::string_view mode    = "neutrino",
      double theta_12          = 0.5836381018669037,
      double theta_23          = 0.8587019919812102,
      double theta_13          = 0.14957471689591406,
      double delta_cp          = 4.084070449666731,
      double dm2_21            = 7.42e-5,
      double dm2_31            = 2.51e-3,
      std::string_view model   = "PREM42");
 
  // ---- Parameter update setters ----
  inline void update_dcp  (double delta_cp)  { chic->update_dcp(delta_cp);     need_update = true; }
  inline void update_th23 (double theta_23)  { chic->update_th23(theta_23);    need_update = true; }
  inline void update_th12 (double theta_12)  { chic->update_th12(theta_12);    need_update = true; }
  inline void update_th13 (double theta_13)  { chic->update_th13(theta_13);    need_update = true; }
  inline void update_dm221(double dm_2_21)   { chic->update_dm221(dm_2_21);    need_update = true; }
  inline void update_dm231(double dm_2_31)   { chic->update_dm231(dm_2_31);    need_update = true; }
  inline void update_density(double rho)     { chic->update_density(rho);      need_update = true; }
  inline void update_Ye(double Ye)           { chic->update_Ye(Ye);            need_update = true; }
 
  // Probability matrix
  Eigen::Matrix3d compute_oscillations(double E, double cos_zenith);
 
  // Return amplitude
  [[nodiscard]] const Eigen::Matrix3cd& get_amplitude() const noexcept { return J; }

private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
 
  // ---- Per-layer cached quantities (energy-dependent) ----
  struct Layer {
    Eigen::Matrix3cd Hs;            // traceless H 
    Eigen::Matrix3cd Hs2;           // traceless H²
    Eigen::Vector3d lambdas;       // eigenvalues
    Eigen::Vector3d diff_lambdas;  // eigenvalue differences
    Eigen::Vector3d prod_lambdas;  // pairwise products of eigenvalues
  };
 
  // ---- Private helpers ----
  void _compute_hamiltonians();                       // fills Layers[]
  void _build_track();                                // fills tracks[]
  void _amplitude();                                  // accumulates J over all layers
  void _layer_amplitude(double L, const Layer& lay);  // single-layer amplitude
 
  // ---- Data members ----
  std::unique_ptr<CHIC>        chic;
  const EarthModel*            Earth = nullptr;       // non-owning, lifetime managed externally
 
  std::vector<Layer>  Layers;   // size = Earth->Nlayers, allocated in constructor
  std::vector<double> tracks, radii2;   // chord half-lengths per layer boundary
 
  // Working matrices — pre-allocated to avoid per-call heap traffic
  Eigen::Matrix3cd J;           // accumulated amplitude
  Eigen::Matrix3cd J_layer;     // single-layer amplitude
  Eigen::Vector3cd exp_eigenvals;
 
  bool   need_update  = true;
  int    deepest      = 0;
  double cos_zenith0, E0;
  std::complex<double> iL;
};
 
#endif // CHICEARTH_H