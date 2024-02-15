#include "Pump.h"

namespace arch {

Pump::Pump() {}

void Pump::setMixture(int mixtureId) {
    this->mixtureId = mixtureId;
}

int Pump::getMixtureId() {
    return mixtureId;
}

}  // namespace arch