#include "Simulator.h"

#include <iostream>
#include <memory>

#include "Channel.h"
#include "Chip.h"
#include "Results.h"
#include "SimulatorImpl.h"

namespace droplet {

Simulator::Simulator() : impl(std::make_unique<SimulatorImpl>()) {}

Simulator::~Simulator() = default;

int Simulator::addChannel(int node0Id, int node1Id, double height, double width, double length) {
    return impl->chip.addChannel(node0Id, node1Id, height, width, length, arch::ChannelType::NORMAL);
}

int Simulator::addMembraneToChannel(int channelId, double height, double width, double poreRadius, double porosity) {
    return impl->chip.addMembraneToChannel(channelId, height, width, poreRadius, porosity);
}

int Simulator::addOrganToMembrane(int membraneId, double height, double width) {
    return impl->chip.addOrganToMembrane(membraneId, height, width);
}

int Simulator::addBypassChannel(int node0Id, int node1Id, double height, double width, double length) {
    return impl->chip.addChannel(node0Id, node1Id, height, width, length, arch::ChannelType::BYPASS);
}

int Simulator::addFlowRatePump(int node0Id, int node1Id, double flowRate) {
    return impl->chip.addFlowRatePump(node0Id, node1Id, flowRate);
}

int Simulator::addPressurePump(int node0Id, int node1Id, double pressure) {
    return impl->chip.addPressurePump(node0Id, node1Id, pressure);
}

void Simulator::addSink(int nodeId) {
    impl->chip.addSink(nodeId);
}

void Simulator::addGround(int nodeId) {
    impl->chip.addGround(nodeId);
}

bool Simulator::checkChipValidity() {
    return impl->chip.isNetworkValid();
}

int Simulator::addFluid(double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient, double saturation) {
    return impl->simulation.addFluid(viscosity, density, concentration, molecularSize, diffusionCoefficient, saturation)->getId();
}

int Simulator::addMixture(std::unordered_map<int, double> fluids) {
    return impl->simulation.addMixture(std::move(fluids));
}

void Simulator::setContinuousPhase(int fluidId) {
    impl->simulation.setContinuousPhase(fluidId);
}

void Simulator::setChangeInputFluid(int fluidId, int pumpId, double injectionTime) {
    impl->simulation.setChangeInputFluid(fluidId, pumpId, injectionTime);
}

void Simulator::setChangeInputMixture(int mixtureId, int pumpId, double injectionTime) {
    impl->simulation.setChangeInputMixture(mixtureId, pumpId, injectionTime);
}

void Simulator::setMaximalAdaptiveTimeStep(double timeStep) {
    impl->simulation.setMaximalAdaptiveTimeStep(timeStep);
}

int Simulator::addDroplet(int fluidId, double volume, double injectionTime, int channelId, double relInjectionPosition) {
    int dropletId = impl->simulation.addDroplet(fluidId, volume)->getId();
    impl->simulation.addInjection(dropletId, injectionTime, channelId, relInjectionPosition);
    return dropletId;
}

void Simulator::setSimulationDuration(double duration) {
    impl->simulation.setSimulationDuration(duration);
}

void Simulator::setSimulationResultTimeStep(double timeStep) {
    impl->simulation.setSimulationResultTimeStep(timeStep);
}

void Simulator::setSimulationCalculationTimeStep(double timeStep) {
    impl->simulation.setSimulationCalculationTimeStep(timeStep);
}

SimulationResult Simulator::simulate() {
    return impl->simulation.simulate();
}

SimulatorImpl::SimulatorImpl() {
    simulation.setChip(&chip);
}

}  // namespace droplet