#include "ChangeInputMixtureEvent.h"

#include "DropletBoundary.h"

namespace sim {

ChangeInputMixtureEvent::ChangeInputMixtureEvent(double time, ChangeInputMixture& injection) : Event(time, 1), injection(injection) {}

void ChangeInputMixtureEvent::performEvent() {
    injection.getInjectionPump()->setMixture(injection.getMixtureId());
    injection.setPerformed(true);
}

}  // namespace sim