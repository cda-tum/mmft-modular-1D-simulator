/**
 * @file IResistanceModel.h
 */

#pragma once

#include "Channel.h"

namespace sim {

class Droplet;

/**
 * @brief Class that specifies the interface of an resistance model. The simulation flow calculations are conducted on the basis of a resistance model.
 */
class IResistanceModel {
  public:
    /**
     * @brief Virtual constructor of a resistance model.
     */
    virtual ~IResistanceModel() {}

    /**
     * @brief Retrieve resistance of the channel itself.
     * @param[in] channel Pointer to channel for which the resistance should be calculated.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel for continuous fluid simulations in which not the continuous phase viscosity is used for calculations. Default -1.0 for cases where it is not used.
     * @return The resistance in the channel itself in Pas/m^3.
     */
    virtual double getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel = -1.0) const = 0;

    /**
     * @brief Retrieve resistance caused by the droplet within the channel.
     * @param[in] channel Pointer to channel in which the droplet currently is.
     * @param[in] droplet Pointer to droplet which causes resistance.
     * @param[in] volumeInsideChannel Volume inside the channel in m^3.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel for continuous fluid simulations in which not the continuous phase viscosity is used for calculations. Default -1.0 for cases where it is not used.
     * @return Resistance caused by the droplet in this channel in Pas/m^3.
     */
    virtual double getDropletResistance(arch::Channel const* const channel, Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel = -1.0) const = 0;
};

}  // namespace sim