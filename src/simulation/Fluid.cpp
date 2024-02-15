#include "Fluid.h"

#include <string>
#include <utility>
#include "Edge.h"

namespace sim {

Fluid::Fluid(int id, double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient, double saturation, int mixtureId) : id(id), viscosity(viscosity), density(density), concentration(concentration), molecularSize(molecularSize), diffusionCoefficient(diffusionCoefficient), saturation(saturation), mixtureId(mixtureId) {}

void Fluid::setName(std::string name) {
    this->name = std::move(name);
}

double Fluid::getSaturation() const {
    return saturation;
}

int Fluid::getId() const {
    return id;
}

std::string Fluid::getName() const {
    return name;
}

double Fluid::getViscosity() const {
    return viscosity;
}

double Fluid::getDensity() const {
    return density;
}

double Fluid::getConcentration() const {
    return concentration;
}

double Fluid::getMolecularSize() const {
    return molecularSize;
}

double Fluid::getDiffusionCoefficient() const {
    return diffusionCoefficient;
}

int Fluid::getMixtureId() const {
    return mixtureId;
}

void Fluid::addMixedFluid(Fluid* fluid) {
    mixedFluids.push_back(fluid);
}

const std::vector<Fluid*>& Fluid::getMixedFluids() const {
    return mixedFluids;
}

}  // namespace sim