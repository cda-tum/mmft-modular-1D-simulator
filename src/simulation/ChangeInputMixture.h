/**
 * @file Injection.h
 */

#pragma once

#include <memory>
#include <string>

#include "ChannelPosition.h"
#include "Droplet.h"
#include "Mixture.h"

namespace sim {

/**
 * @brief Class that contains all paramaters necessary to conduct an injection.
 */
class ChangeInputMixture {
  private:
    const int id;               ///< Unique identifier of an injection.
    std::string name = "";      ///< Name of the injection.
    double injectionTime;       ///< Time at which the injection should start in s elapsed since the start of the simulation.
    int mixtureId;              //< Id of mixture to be injected.
    arch::Pump* injectionPump;  ///< Pump at which the fluid is injected.
    bool performed = false;     ///< Information if the change of the input mixture was already performed or not.

  public:
    /**
     * @brief Create an injection.
     * @param[in] id Unique identifier of an injection.
     * @param[in] mixture Pointer to the mixture to be injected.
     * @param[in] injectionPump Pump from which the fluid is injection
     * @param[in] injectionTime Time at which the fluid injection should start, in s elapsed since the start of the simulation in s.
     */
    ChangeInputMixture(int id, int mixtureId, arch::Pump* injectionPump, double injectionTime);

    /**
     * @brief Set name of injection.
     * @param[in] name Name of injection.
     */
    void setName(std::string name);

    /**
     * @brief Retrieve unique identifier of injection.
     * @return Unique identifier of injection.
     */
    int getId() const;

    /**
     * @brief Retireve name of injection.
     * @return Name of injection.
     */
    std::string getName() const;

    /**
     * @brief Retrieve time at which the droplet should be injected.
     * @return Time in s elapsed since the start of the simulation at which the injection of the fluid should take place in s.
     */
    double getInjectionTime() const;

    /**
     * @brief Retrieve pump from which the injection should take place.
     * @return Pump from which the injection should take place. 
     */
    arch::Pump* getInjectionPump();

    /**
     * @brief Retrieve pointer to the fluid that should be injected.
     * @return Mixture to be injected.
     */
    int getMixtureId();

    /**
     * @brief Update the state if the input mixture change was already performed or not.
     * 
     * @param performed If the input mixture change was performed or not.
     */
    void setPerformed(bool performed);

    /**
     * @brief Get information if the change of the input mixture was already performed or not.
     * 
     * @return true: Change was already performed.
     * @return false: Change was not yet performed. 
     */
    bool wasPerformed();
};

}  // namespace sim