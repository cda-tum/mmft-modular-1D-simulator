/**
 * @file ResistanceModels.h
 */

#pragma once

#include "Fluid.h"
#include "IResistanceModel.h"

namespace sim {

/**
 * @brief Class that defines the functionality of the 1D resistance model.
 */
class ResistanceModel0 : public sim::IResistanceModel {
  private:
    double continuousPhaseViscosity;  ///< The viscosity of the continuous phase in Pas.

  public:
    /**
     * @brief Instantiate the resistance model.
     * @param[in] continuousPhaseViscosity The viscosity of the continuous phase in Pas.
     */
    ResistanceModel0(double continuousPhaseViscosity);

    /**
     * @brief Calculate and returns the resistance of the channel itself.
     * @param[in] channel A pointer to the channel for which the resistance should be calculated.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Not used in this resistance model and therefore default -1.0.
     * @return The resistance of the channel itself in Pas/m^3.
     */
    double getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const override;

    /**
     * @brief Compute the a factor.
     * @param[in] width Width of the channel in m.
     * @param[in] height Height of the channel in m.
     * @return The a factor.
     */
    double computeFactorA(double width, double height) const;

    /**
     * @brief Retrieve the resistance caused by the specified droplet in the specified channel.
     * @param[in] channel Pointer to channel for which the droplet resistance should be calculated.
     * @param[in] droplet Pointer to droplet that causes the droplet resistance in the channel.
     * @param[in] volumeInsideChannel The volume inside the channel in m^3.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Not used in this resistance model and therefore default -1.0.
     * @return Resistance caused by the droplet in the channel in Pas/m^3.
     */
    double getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const override;
};

/**
 * @brief Class that defines a test resistance model.
 */
class ResistanceModel1 : public sim::IResistanceModel {
  public:
    /**
     * @brief Instantiate the resistance model.
     */
    ResistanceModel1();

    /**
     * @brief Calculate and returns the resistance of the channel itself.
     * @param[in] channel Channel for which the resistance should be calculated.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Not used in this resistance model and therefore default -1.0.
     * @return The resistance of the channel in Pas/m^3.
     */
    double getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const override;

    /**
     * @brief Retrieve the resistance caused by the specified droplet in the specified channel.
     * @param[in] channel Pointer to channel for which the droplet resistance should be calculated.
     * @param[in] droplet Pointer to droplet that causes the droplet resistance in the channel.
     * @param[in] volumeInsideChannel The volume inside the channel in m^3.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Not used in this resistance model and therefore default -1.0.
     * @return Resistance caused by the droplet in the channel in Pas/m^3.
     */
    double getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const override;
};

/**
 * @brief Class that defines the functionality of the 1D resistance model for continuous fluid simulation.
 */
class ResistanceModel2 : public sim::IResistanceModel {
  private:
    double continuousPhaseViscosity;  ///< The average viscosity of continuous fluids in a channel in Pas.

  public:
    /**
     * @brief Instantiate the resistance model.
     * @param[in] continuousPhaseViscosity ///< The viscosity of the continuous phase in Pas.
     */
    ResistanceModel2(double continuousPhaseViscosity);

    /**
     * @brief Calculate and returns the resistance of the channel itself.
     * @param[in] channel A pointer to the channel for which the resistance should be calculated.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Used instead of viscosity of continuous phase.
     * @return The resistance of the channel itself in Pas/m^3.
     */
    double getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const override;

    /**
     * @brief Compute the a factor.
     * @param[in] width Width of the channel in m.
     * @param[in] height Height of the channel in m.
     * @return The a factor.
     */
    double computeFactorA(double width, double height) const;

    /**
     * @brief Retrieve the resistance caused by the specified droplet in the specified channel.
     * @param[in] channel Pointer to channel for which the droplet resistance should be calculated.
     * @param[in] droplet Pointer to droplet that causes the droplet resistance in the channel.
     * @param[in] volumeInsideChannel The volume inside the channel in m^3.
     * @param[in] avgViscosityInChannel Avergage viscosity of continuous fluids in channel. Used instead of continuous phase viscosity.
     * @return Resistance caused by the droplet in the channel in Pas/m^3.
     */
    double getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const override;
};

}  // namespace sim