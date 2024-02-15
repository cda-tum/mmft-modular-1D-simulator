#include "StopSimulationEvent.h"

namespace sim {

StopSimulationEvent::StopSimulationEvent(double time, Simulation& simulation) : Event(time, 2), simulation(simulation) {}

void StopSimulationEvent::performEvent() {
    // this event does nothing and only ensures a minimal time step, so a new flow state is computed again
    simulation.stopSimulation();
}

}  // namespace sim