#include "MembraneResistanceModels.h"

#include <math.h>
#include <cmath>
#include "Fluid.h"
#include "Mixture.h"

namespace sim {

// ### MembraneResistanceModel0 ###
// Ishahak2020
MembraneResistanceModel0::MembraneResistanceModel0() {}

double MembraneResistanceModel0::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return getPoreResistance(membrane, fluid) / (area * membrane->getPorosity());
    //return 1 / getPermeabilityParameter(membrane);
}

double MembraneResistanceModel0::getPermeabilityParameter(const arch::Membrane* const membrane) const {
    return (M_PI * membrane->getPorosity() * pow(membrane->getPoreRadius(), 4)) / 8;
}

double MembraneResistanceModel0::getPoreResistance(arch::Membrane const* const membrane, Fluid const* const fluid) const {
    return (8 * fluid->getViscosity() * membrane->getHeight()) / (M_PI * pow(membrane->getPoreRadius(), 4));
}

// ### MembraneResistanceModel1 ###
// VanDersari2011
MembraneResistanceModel1::MembraneResistanceModel1() {}

double MembraneResistanceModel1::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    double deviceAdjustment = 1.9;
    return deviceAdjustment * getPoreDensityDependentResistance(membrane, area, fluid->getDiffusionCoefficient()) + (getPoreResistance(membrane, fluid->getDiffusionCoefficient()) + getPoreExitResistance(membrane, fluid->getDiffusionCoefficient())) / membrane->getNumberOfPores(area);
}

double MembraneResistanceModel1::getPoreExitResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return membrane->getHeight() / (diffusionCoefficient * M_PI * pow(membrane->getPoreRadius(), 2));
}

double MembraneResistanceModel1::getPoreResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (4 * diffusionCoefficient * membrane->getPoreRadius());
}

double MembraneResistanceModel1::getPoreDensityDependentResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const {
    double proportionalityConstant = 5.3;
    return proportionalityConstant / (diffusionCoefficient * area / membrane->getWidth() * M_PI);
}

// ### MembraneResistanceModel2 ###
//Snyder2011
MembraneResistanceModel2::MembraneResistanceModel2() {}

double MembraneResistanceModel2::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    double moleculeToPoreRadius = (fluid->getMolecularSize() / 2) / (membrane->getPoreRadius());
    double reducedDiffusionCoefficient = (-2.81903 * moleculeToPoreRadius + 0.270788 * moleculeToPoreRadius + 1.1015 * moleculeToPoreRadius - 0.435933 * moleculeToPoreRadius) * fluid->getDiffusionCoefficient();
    return 1 / (1 / getPoreDiscoveryResistance(membrane, area, fluid->getDiffusionCoefficient()) + 1 / getTransmembraneResistance(membrane, reducedDiffusionCoefficient, area));
}

double MembraneResistanceModel2::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double area, double freeDiffusionCoefficient) const {
    return (2 * freeDiffusionCoefficient / membrane->getArea()) * membrane->getPoreRadius() * membrane->getNumberOfPores(area);
}

double MembraneResistanceModel2::getTransmembraneResistance(arch::Membrane const* const membrane, double reducedDiffusionCoefficient, double area) const {
    return 1 / (area * membrane->getHeight()) * membrane->getNumberOfPores(area) * reducedDiffusionCoefficient * M_PI * pow(membrane->getPoreRadius(), 2);
}

// ### MembraneResistanceModel3 ###
//Berg1993, assume sphere
MembraneResistanceModel3::MembraneResistanceModel3() {}

double MembraneResistanceModel3::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return getPorePassageResistance(membrane, fluid->getDiffusionCoefficient(), area) + getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) / membrane->getNumberOfPores(area);
}

double MembraneResistanceModel3::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (4 * membrane->getPoreRadius() * diffusionCoefficient);
}

double MembraneResistanceModel3::getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient, double area) const {
    return 1 / (4 * M_PI * (std::sqrt(area / M_PI)) * diffusionCoefficient);
}

// ### MembraneResistanceModel4 ###
//Chung2018
MembraneResistanceModel4::MembraneResistanceModel4() {}

double MembraneResistanceModel4::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return (getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) + getPorePassageResistance(membrane, fluid->getDiffusionCoefficient())) / membrane->getNumberOfPores(area);
}

double MembraneResistanceModel4::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (4 * membrane->getPoreRadius() * diffusionCoefficient);
}

double MembraneResistanceModel4::getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return membrane->getHeight() / (4 * M_PI * std::pow(membrane->getPoreRadius(), 2) * diffusionCoefficient);
}

// ### MembraneResistanceModel5 ###
//Chung2018 Appendix 1 Calculations
MembraneResistanceModel5::MembraneResistanceModel5() {}

double MembraneResistanceModel5::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return (getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) + getPorePassageResistance(membrane, fluid->getDiffusionCoefficient())) / membrane->getNumberOfPores(area);
}

double MembraneResistanceModel5::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (4 * membrane->getPoreRadius() * 2 * diffusionCoefficient);
}

double MembraneResistanceModel5::getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return membrane->getHeight() / (4 * M_PI * std::pow(membrane->getPoreRadius() * 2, 2) * diffusionCoefficient);
}

// ### MembraneResistanceModel6 ###
//Ronaldson-Bouchard 2022
MembraneResistanceModel6::MembraneResistanceModel6() {}

double MembraneResistanceModel6::getMembraneResistance(arch::Membrane const* const membrane, sim::Fluid const* const fluid, double area) const {
    return (0.5 * membrane->getChannel()->getHeight() / fluid->getDiffusionCoefficient() + membrane->getHeight() / diffusionFactorMembrane + membrane->getOrgan()->getHeight() * 0.5 / fluid->getDiffusionCoefficient());
}

// ### MembraneResistanceModel7 ###
//Berg1993, permeability increases twice as fast for small number of N
MembraneResistanceModel7::MembraneResistanceModel7() {}

double MembraneResistanceModel7::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return 0.5 * (getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) / membrane->getNumberOfPores(area) + getPorePassageResistance(membrane, area, fluid->getDiffusionCoefficient()));
}

double MembraneResistanceModel7::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (2 * membrane->getPoreRadius() * diffusionCoefficient);
}

double MembraneResistanceModel7::getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const {
    return membrane->getHeight() / (diffusionCoefficient * area);
}

// ### MembraneResistanceModel8 ###
//Berg1993, planar barrier
MembraneResistanceModel8::MembraneResistanceModel8() {}

double MembraneResistanceModel8::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return (getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) / membrane->getNumberOfPores(area) + getPorePassageResistance(membrane, area, fluid->getDiffusionCoefficient()));
}

double MembraneResistanceModel8::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (2 * membrane->getPoreRadius() * diffusionCoefficient);
}

double MembraneResistanceModel8::getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const {
    return membrane->getHeight() / (diffusionCoefficient * area);
}

// ### MembraneResistanceModel9 ###
//Berg1993, planar barrier with disc-absorber
MembraneResistanceModel9::MembraneResistanceModel9() {}

double MembraneResistanceModel9::getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const {
    return (getPoreDiscoveryResistance(membrane, fluid->getDiffusionCoefficient()) / membrane->getNumberOfPores(area) + getPorePassageResistance(membrane, area, fluid->getDiffusionCoefficient()));
}

double MembraneResistanceModel9::getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const {
    return 1 / (4 * membrane->getPoreRadius() * diffusionCoefficient);
}

double MembraneResistanceModel9::getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const {
    return membrane->getHeight() / (diffusionCoefficient * area);
}

}  // namespace sim