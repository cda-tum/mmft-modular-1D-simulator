/**
 * @file IMembraneResistanceModel.h
 */

#pragma once

#include "Membrane.h"
#include "Mixture.h"

namespace sim {

/**
 * @brief Class that specifies the interface of an membrane resistance model. The membrane resistance calculations are conducted on the basis of a resistance model.
 */
class IMembraneResistanceModel {
  public:
    /**
     * @brief Virtual constructor of a resistance model.
     */
    virtual ~IMembraneResistanceModel() {}

    /**
     * @brief Get the resistance caused by the membrane based on the specific resistance model.
     * 
     * @param membrane Pointer to the membrane for which the resistance should be calculated.
     * @param fluid Species for which the resistance should be calculated.
     * @return The resistance of this membrane to this fluid.
     */
    virtual double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const = 0;
};

}  // namespace sim