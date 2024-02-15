/**
 * @file ChangeInputMixtureEvent.h
 */

#pragma once

#include "ChangeInputMixture.h"
#include "Event.h"

namespace sim {

/**
 * @brief Class for an injection event that takes place when a continuous fluid injection is started.
 */
class ChangeInputMixtureEvent : public Event {
  private:
    ChangeInputMixture& injection;  ///< Specifies if the injection event.

  public:
    /**
     * @brief Definies an injection event to take place at a certain time.
     * @param[in] time The time at which the event should take place in s elapsed since the start of the simulation.
     * @param[in,out] injection A class containing all details necessary to conduct an continuous injection event.
     */
    ChangeInputMixtureEvent(double time, ChangeInputMixture& injection);

    /**
     * @brief Conducts the injection event.
     */
    void performEvent() override;
};

}  // namespace sim