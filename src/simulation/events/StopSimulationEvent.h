/**
 * @file StopSimulationEvent.h
 */
#pragma once

#include "Event.h"
#include "Simulation.h"

namespace sim {

/**
 * @brief Class to trigger the calculation of the simulation parameters after a minimal time step.
 */
class StopSimulationEvent : public Event {
  private:
    Simulation& simulation;  ///< Simulation class
  public:
    /**
     * @brief Construct class to schedule a minimal tim estep event.
     * @param[in] time Time after minimal time step passed in s elapsed since the start of the simulation.
     */
    StopSimulationEvent(double time, Simulation& simulation);

    /**
     * @brief Do nothing except for logging the event. As the event exists, the simulation will be forwarded to this time point in the simulation algorithm and therefore it is ensured that the simulation parameters at this point in time are calculated.
     */
    void performEvent() override;
};

}  // namespace sim