#include "ChangeInputMixture.h"

#include <memory>
#include <utility>

#include "ChannelPosition.h"

namespace sim {

ChangeInputMixture::ChangeInputMixture(int id, int mixtureId, arch::Pump* injectionPump, double injectionTime) : id(id), mixtureId(mixtureId), injectionPump(injectionPump), injectionTime(injectionTime) {}

void ChangeInputMixture::setName(std::string name) {
    name = std::move(name);
}

int ChangeInputMixture::getId() const {
    return id;
}

std::string ChangeInputMixture::getName() const {
    return name;
}

double ChangeInputMixture::getInjectionTime() const {
    return injectionTime;
}

arch::Pump* ChangeInputMixture::getInjectionPump() {
    return injectionPump;
}

int ChangeInputMixture::getMixtureId() {
    return mixtureId;
}

void ChangeInputMixture::setPerformed(bool performed) {
    this->performed = performed;
}

bool ChangeInputMixture::wasPerformed() {
    return performed;
}

}  // namespace sim