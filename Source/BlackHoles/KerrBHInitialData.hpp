/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KERRBHINITIALDATA_HPP_
#define KERRBHINITIALDATA_HPP_

#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include <array>
#include <cmath>

//! Class which computes the Kerr initial conditions per arXiv 1401.1548
class KerrBHInitialData
{
 
   public:
    //! Stuct for the params of the Kerr BHInitialData
    struct params_t
    {
        amrex::Real mass; //!<< The mass of the Kerr BH
        std::array<amrex::Real, AMREX_SPACEDIM>
            center;       //!< The center of the Kerr BH
        amrex::Real spin; //!< The spin param a = J/M, so 0 <= |a| <= M
        std::array<amrex::Real, AMREX_SPACEDIM>
            spin_direction; //!< The direction of the spin of the Kerr BH
    };

  protected:
    amrex::Real m_dx;
    params_t m_params;
    Tensor::Rank2 m_R;

  public:
    KerrBHInitialData(params_t a_params, amrex::Real a_dx)
        : m_dx(a_dx), m_params(a_params)

    {
        // check this spin param is sensible
        if (std::abs(m_params.spin) > m_params.mass)
        {
            amrex::Abort("The spin parameter must satisfy |a| <= M");
        }
        Tensor::Rank1 spin_dir{
            m_params.spin_direction[0],
            m_params.spin_direction[1], 
            m_params.spin_direction[2]};

        amrex::Real spin_dir_norm = sqrt(spin_dir(0) * spin_dir(0) + spin_dir(1) * spin_dir(1) +
                                      spin_dir(2) * spin_dir(2));
        if (spin_dir_norm <1e-14)
        {
            amrex::Abort("The spin direction must be non-zero");
        }
        FOR(i)
        {
            spin_dir(i) /= spin_dir_norm;
        }

        Tensor::Rank1 z_dir = {0.0, 0.0, 1.0};

        m_R = CoordinateTransformations::rotation_matrix(spin_dir, z_dir);
    }

    void operator()(int ix, int iy, int iz,
           const amrex::Array4<amrex::Real> &state) const;
  protected:
    //! Function which computes the components of the metric in spherical coords
    AMREX_FORCE_INLINE 
    void compute_kerr(
    Tensor::Sym12Rank2 &spherical_g,
    Tensor::Sym12Rank2 &spherical_K,
    Tensor::Rank1 &spherical_shift,
    amrex::Real &kerr_lapse,
    const Coordinates coords) const;
};

#include "KerrBHInitialData.impl.hpp"

#endif /* KERRBHINITIALDATA_HPP_ */
