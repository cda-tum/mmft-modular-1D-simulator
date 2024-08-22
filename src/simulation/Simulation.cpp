#include "Simulation.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "BoundaryHeadEvent.h"
#include "BoundaryTailEvent.h"
#include "ChangeInputMixture.h"
#include "ChangeInputMixtureEvent.h"
#include "Channel.h"
#include "ChannelPosition.h"
#include "Droplet.h"
#include "DropletBoundary.h"
#include "Fluid.h"
#include "IFlowRatePump.h"
#include "INode.h"
#include "IPressurePump.h"
#include "IResistance.h"
#include "Injection.h"
#include "InjectionEvent.h"
#include "MergeBifurcationEvent.h"
#include "MergeChannelEvent.h"
#include "NodalAnalysis.h"
#include "ResistanceModels.h"
#include "Results.h"
#include "StopSimulationEvent.h"
#include "TimeStepEvent.h"

#include <cassert>

namespace sim {

Simulation::Simulation() {}

void Simulation::setChip(arch::Chip* chip) {
    this->chip = chip;
}

arch::Chip* Simulation::getChip() {
    return chip;
}

void Simulation::setContinuousPhase(int fluidId) {
    continuousPhaseMixtureId = fluids.at(fluidId)->getMixtureId();
}

int Simulation::getMixtureOfFluid(int fluidId) {
    for (auto& mixture : mixtures) {
        if (mixture.getConcentrationOfFluid(fluidId) == 1.0) {
            return mixture.getId();
        }
    }
    throw std::invalid_argument("No mixture with 1.0 concentration of this fluid exists.");  //should never happen as mixture is added when fluid is added.
}

void Simulation::fillChipWithContinuousPhaseMixture() {
    for (auto& [channelId, channel] : chip->getChannels()) {
        std::deque<std::pair<int, double>> mixture_pos = {std::make_pair(continuousPhaseMixtureId, 1.0)};
        mixturesInEdge.try_emplace(channelId, mixture_pos);
    }

    for (auto& [membranId, membrane] : chip->getMembranes()) {
        std::deque<std::pair<int, double>> mixture_pos = {std::make_pair(continuousPhaseMixtureId, 1.0)};
        mixturesInEdge.try_emplace(membranId, mixture_pos);
        mixturesInEdge.try_emplace(membrane->getOrgan()->getId(), mixture_pos);
    }
    for (auto& [pumpId, pump] : chip->getFlowRatePumps()) {
        pump->setMixture(continuousPhaseMixtureId);
        std::deque<std::pair<int, double>> mixture_pos = {std::make_pair(continuousPhaseMixtureId, 1.0)};
        mixturesInEdge.try_emplace(pumpId, mixture_pos);
    }
    for (auto& [pumpId, pump] : chip->getPressurePumps()) {
        pump->setMixture(continuousPhaseMixtureId);
        std::deque<std::pair<int, double>> mixture_pos = {std::make_pair(continuousPhaseMixtureId, 1.0)};
        mixturesInEdge.try_emplace(pumpId, mixture_pos);
    }
}

ChangeInputMixture* Simulation::setChangeInputFluid(int fluidId, int pumpId, double injectionTime) {
    auto fluid = getFluid(fluidId);
    auto mixture = fluid->getMixtureId();
    auto injectionPump = chip->getPump(pumpId);

    auto id = changeInputMixtures.size();

    auto result = changeInputMixtures.insert_or_assign(id, std::make_unique<ChangeInputMixture>(id, mixture, injectionPump, injectionTime));

    return result.first->second.get();
}

ChangeInputMixture* Simulation::setChangeInputMixture(int mixtureId, int pumpId, double injectionTime) {
    auto injectionPump = chip->getPump(pumpId);

    auto id = changeInputMixtures.size();

    auto result = changeInputMixtures.insert_or_assign(id, std::make_unique<ChangeInputMixture>(id, mixtureId, injectionPump, injectionTime));

    return result.first->second.get();
}

ChangeInputMixture* Simulation::getChangeInputFluid(int injectionId) {
    return changeInputMixtures.at(injectionId).get();
}

void Simulation::setResistanceModel(ResistanceModel modelName) {
    this->resistanceModelName = modelName;
}

void Simulation::setMembraneResistanceModel(MembraneResistanceModel modelName) {
    this->membraneResistanceModelName = modelName;
}

void Simulation::setMaximalAdaptiveTimeStep(double timeStep) {
    maximalAdaptiveTimeStep = timeStep;
}

Fluid* Simulation::addFluid(double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient, double saturation) {
    int id = fluids.size();
    int mixtureId = mixtures.size();

    auto result = fluids.insert_or_assign(id, std::make_unique<Fluid>(id, viscosity, density, concentration, molecularSize, diffusionCoefficient, saturation, mixtureId));

    addMixture(std::unordered_map<int, double>{{id, concentration}}, mixtureId);

    return result.first->second.get();
}

Fluid* Simulation::getFluid(int fluidId) {
    return fluids.at(fluidId).get();
}

Mixture* Simulation::getMixture(int mixtureId) {
    return &mixtures.at(mixtureId);
}

int Simulation::addMixture(std::unordered_map<int, double>&& fluidConcentrations, int id) {
    if (id == -1) {  // if default id parameter is used
        id = mixtures.size();
    }

    double viscosity = 0;
    double density = 0;
    double largestMolecularSize = 0;
    for (auto& [fluidId, concentration] : fluidConcentrations) {
        viscosity += fluids.at(fluidId)->getViscosity() * concentration;
        density += fluids.at(fluidId)->getDensity() * concentration;
        if (fluids.at(fluidId)->getMolecularSize() > largestMolecularSize) {
            largestMolecularSize = fluids.at(fluidId)->getMolecularSize();
        }
    }

    mixtures.emplace_back(id, std::move(fluidConcentrations), viscosity, density, largestMolecularSize);

    return id;
}

Droplet* Simulation::addDroplet(int fluidId, double volume) {
    auto id = droplets.size();
    auto fluid = fluids.at(fluidId).get();

    auto result = droplets.insert_or_assign(id, std::make_unique<Droplet>(id, volume, fluid));

    return result.first->second.get();
}

void Simulation::setSimulationDuration(double duration) {
    this->simulationDuration = duration;
}

void Simulation::setSimulationResultTimeStep(double timeStep) {
    this->simulationResultTimeStep = timeStep;
}

void Simulation::setSimulationCalculationTimeStep(double timeStep) {
    this->minSimulationCalcTimeStep = timeStep;
}

Droplet* Simulation::getDroplet(int dropletId) {
    return droplets.at(dropletId).get();
}

Droplet* Simulation::getDropletAtNode(int nodeId) {
    // loop through all droplets
    for (auto& [id, droplet] : droplets) {
        // do not consider droplets which are not inside the network
        if (droplet->getDropletState() != DropletState::NETWORK) {
            continue;
        }

        // do not consider droplets inside a single channel (because they cannot span over a node)
        if (droplet->isInsideSingleChannel()) {
            continue;
        }

        // check if a boundary is connected with the reference node and if yes then return the droplet immediately
        auto connectedBoundaries = droplet->getConnectedBoundaries(nodeId);
        if (connectedBoundaries.size() != 0) {
            return droplet.get();
        }

        // check if a fully occupied channel is connected with the reference node and if yes then return the droplet immediately
        auto connectedFullyOccupiedChannels = droplet->getConnectedFullyOccupiedChannels(nodeId);
        if (connectedFullyOccupiedChannels.size() != 0) {
            return droplet.get();
        }
    }

    // if nothing was found than return nullptr
    return nullptr;
}

Injection* Simulation::addInjection(int dropletId, double injectionTime, int channelId, double injectionPosition) {
    auto id = injections.size();
    auto droplet = droplets.at(dropletId).get();
    auto channel = chip->getChannel(channelId);

    // --- check if injection is valid ---
    // for the injection the head and tail of the droplet must lie inside the channel (the volume of the droplet must be small enough)
    // the droplet length is a relative value between 0 and 1
    double dropletLength = droplet->getVolume() / channel->getVolume();
    // channel must be able to fully contain the droplet
    if (dropletLength >= 1.0) {
        throw std::invalid_argument("Injection of droplet " + droplet->getName() + " into channel " + channel->getName() + " is not valid. Channel must be able to fully contain the droplet.");
    }
    // compute tail and head position of the droplet
    double tail = (injectionPosition - dropletLength / 2);
    double head = (injectionPosition + dropletLength / 2);
    // tail and head must not be outside the channel (can happen when the droplet is injected at the beginning or end of the channel)
    if (tail < 0 || head > 1.0) {
        throw std::invalid_argument("Injection of droplet " + droplet->getName() + " is not valid. Tail and head of the droplet must lie inside the channel " + channel->getName() + ". Consider to set the injection position in the middle of the channel.");
    }

    auto result = injections.insert_or_assign(id, std::make_unique<Injection>(id, droplet, injectionTime, channel, injectionPosition));

    return result.first->second.get();
}

Injection* Simulation::getInjection(int injectionId) {
    return injections.at(injectionId).get();
}

droplet::SimulationResult Simulation::simulate() {
    // ##########
    // Initialize
    // ##########

    droplet::SimulationResult result;
    storeSimulationParameters(result);

    // initialize the simulation
    initialize();

    // get nodes, channels, pumps
    // nodes (<nodeId, <Node, matrixId>>): key is the nodeId and the tuple consists of the node and the correspondng matrixId (the groundNode has matrixId -1)
    // this is done, because the nodal analysis needs a nodeId that represents the node's index in the matrix (hence, the matrixId), which is usually not the same as the nodeId
    // the node itself inside the tuple is not really necessary, but is convenient when writing all pressure values to each node, after solving the matrix
    std::unordered_map<int, std::tuple<nodal::INode*, int>> nodes;
    int matrixId = 0;
    for (auto& [nodeId, node] : chip->getNodes()) {
        // if node is ground node set the matrixId to -1
        if (nodeId == chip->getGroundId()) {
            nodes.insert_or_assign(nodeId, std::make_tuple(node.get(), -1));
        } else {
            nodes.insert_or_assign(nodeId, std::make_tuple(node.get(), matrixId));
            matrixId++;
        }
    }
    std::vector<nodal::IResistance*> channels;
    for (auto& [key, channel] : chip->getChannels()) {
        channels.push_back(channel.get());
    }

    // ##########
    // Simulation Loop
    // ##########
    // * update droplet resistances
    // * conduct nodal analysis
    // * update droplets (flow rates of boundaries)
    // * compute events
    // * search for next event (break if no event is left)
    // * move droplets
    // * perform event
    double simulationResultTimeCounter = 0;
    simulationInProgress = true;
    while (simulationInProgress) {
        // update droplet resistances (in the first iteration no  droplets are inside the network)
        updateChannelResistances();
        updateDropletResistances();

        std::vector<nodal::IFlowRatePump*> flowRatePumps;
        for (auto& [key, pump] : chip->getFlowRatePumps()) {
            flowRatePumps.push_back(pump.get());
        }
        std::vector<nodal::IPressurePump*> pressurePumps;
        for (auto& [key, pump] : chip->getPressurePumps()) {
            pressurePumps.push_back(pump.get());
        }

        // compute nodal analysis
        nodal::conductNodalAnalysis(nodes, channels, pressurePumps, flowRatePumps);

        // update droplets, i.e., their boundary flow rates
        // loop over all droplets
        dropletsAtBifurcation = false;
        for (auto& [key, droplet] : droplets) {
            // only consider droplets inside the network
            if (droplet->getDropletState() != DropletState::NETWORK) {
                continue;
            }

            // set to true if droplet is at bifurcation
            if (droplet->isAtBifurcation()) {
                dropletsAtBifurcation = true;
            }

            // compute the average flow rates of all boundaries, since the inflow does not necessarily have to match the outflow (qInput != qOutput)
            // in order to avoid an unwanted increase/decrease of the droplet volume an average flow rate is computed
            // the actual flow rate of a boundary is then determined accordingly to the ratios of the different flowRates inside the channels
            droplet->updateBoundaries(*chip, slipFactor);
        }

        // compute internal minimal timestep
        double minTime = std::numeric_limits<double>::max();
        for (auto& [key, channel] : chip->getChannels()) {
            auto flowrate = channel->getFlowRate();
            auto length = channel->getLength();
            auto time = length * channel->getArea() / std::abs(flowrate);
            if (time < minTime) {
                minTime = time;
            }
        }
        //std::cout << "minTime " << minTime << std::endl;
        if (minSimulationCalcTimeStep > simulationResultTimeStep) {
            minSimulationCalcTimeStep = simulationResultTimeStep;
        }
        if (minTime < minSimulationCalcTimeStep) {
            internalSimulationTimeStep = minTime - 0.5e-10;  // 0.5e-10 factor to avoid pos larger than 1.0 due to rounding errors
        } else {
            internalSimulationTimeStep = minSimulationCalcTimeStep;
        }
        assert(internalSimulationTimeStep <= simulationResultTimeStep && internalSimulationTimeStep <= minSimulationCalcTimeStep);

        // merging:
        // 1) Two boundaries merge inside a single channel:
        //      1) one boundary is faster than the other and would "overtake" the boundary
        //          * theoretically it would also be possible to introduce a merge distance.
        //            however, then it gets complicated when the actual merging of the droplets happen, since all other boundaries of the droplets need to moved in order to conserve the droplet volume (i.e., compensate for the volume between the droplets due to the merge distance).
        //            this movement of all other boundaries would then may trigger other events => it gets complicated^^
        //      2) the two boundaries flow in different directions:
        //          * this case shouldn't happen in this implementation now, but it might be possible since boundary flow rates are not bound to channel flow rates.
        //            however, currently the boundaries should still have the same direction as the channel flow rates and, thus, the boundaries cannot have opposite flow directions in a single channel
        // 2) A boundary (which currently operates as a droplet head) can merge into another droplet when it switches channels (basically when a BoundaryHeadEvent occurs):
        //      * possible solution is to normally compute a BoundaryHeadEvent and then check if it is actually a merge event (in case of a merge event the channels don't have to be switched since the boundary is merged with the other droplet)

        // store simulation results of current state
        if (simulationResultTimeCounter <= 0) {
            storeSimulationResults(result);
            if (changeInputMixtures.empty() || !injections.empty()) {
                simulationResultTimeCounter = 0;
            } else {  // for continuous fluid simulation without droplets only store simulation time steps
                simulationResultTimeCounter = simulationResultTimeStep;
            }
        }

        // compute events
        auto events = computeEvents(simulationResultTimeCounter);

        // sort events
        // closest events in time with the highest priority come first
        std::sort(events.begin(), events.end(), [](auto& a, auto& b) {
            if (a->getTime() == b->getTime()) {
                return a->getPriority() < b->getPriority();  // ascending order (the lower the priority value, the higher the priority)
            }
            return a->getTime() < b->getTime();  // ascending order
        });

        // get next event or break loop, if no events remain
        Event* nextEvent = nullptr;
        if (events.size() != 0) {
            nextEvent = events[0].get();
        } else {
            break;
        }

        auto nextEventTime = nextEvent->getTime();
        currTime += nextEventTime;
        simulationResultTimeCounter -= nextEventTime;

        if (nextEventTime > 0.0) {
            // move droplets until event is reached
            moveDroplets(nextEventTime);

            // update fluid concentrations in channels
            //updateFluidConcentrationsInChannels(nextEventTime);
            if (!changeInputMixtures.empty()) {
                calculateNewMixtures(nextEventTime);
            }
        }

        // perform event (inject droplet, move droplet to next channel, block channel, merge channels)
        // only one event at a time is performed in order to ensure a correct state
        // hence, it might happen that the next time for an event is 0 (e.g., when  multiple events happen at the same time)
        assert(nextEventTime >= 0.0);
        nextEvent->performEvent();
    }

    storeSimulationResults(result);  //store simulation parameters once more at the end of the simulation (even if this is not part of the simulation result timestep)

    return result;
}

void Simulation::initialize() {
    auto* continuousPhaseMixture = &mixtures.at(continuousPhaseMixtureId);
    // set resistance model
    if (resistanceModelName == ResistanceModel::ONE_D_MODEL) {
        resistanceModel = std::make_unique<ResistanceModel0>(continuousPhaseMixture->getViscosity());
    } else if (resistanceModelName == ResistanceModel::TEST_MODEL) {
        resistanceModel = std::make_unique<ResistanceModel1>();
    } else if (resistanceModelName == ResistanceModel::ONE_D_CONTINUOUS_MODEL) {
        resistanceModel = std::make_unique<ResistanceModel2>(continuousPhaseMixture->getViscosity());
    }

    // set membrane resistance model
    if (membraneResistanceModelName == MembraneResistanceModel::MODEL_0) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel0>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_1) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel1>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_2) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel2>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_3) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel3>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_4) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel4>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_5) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel5>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_6) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel6>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_7) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel7>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_8) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel8>();
    } else if (membraneResistanceModelName == MembraneResistanceModel::MODEL_9) {
        membraneResistanceModel = std::make_unique<MembraneResistanceModel9>();
    }

    // compute channel resistances
    for (auto& [key, channel] : chip->getChannels()) {
        double resistance = resistanceModel->getChannelResistance(channel.get());
        channel->setChannelResistance(resistance);
        channel->setDropletResistance(0.0);
    }

    fillChipWithContinuousPhaseMixture();
}

double Simulation::getAverageViscosityInChannel(int channelId) {
    double avgViscosity = 0.0;
    double prevPos = 0.0;
    for (auto it = mixturesInEdge.at(channelId).crbegin(); it != mixturesInEdge.at(channelId).crend(); it++) {
        auto& [mixtureId, pos] = *it;
        //std::cout << "viscosity" << mixtures.at(mixtureId).getViscosity() << "pos " << pos << "prevPos " << prevPos << "(pos - prevPos)" << pos - prevPos << std::endl;
        avgViscosity += mixtures.at(mixtureId).getViscosity() * (pos - prevPos);
        prevPos = pos;
    }
    //std::cout << "avgViscosity" << avgViscosity << std::endl;
    return avgViscosity;
}

void Simulation::updateChannelResistances() {
    // the channel resistances are only updated if the simulation resistance model considers the average viscosity of all fluids in a channel (instead of just the continuous phase viscosity that does not change during the simulation)
    if (resistanceModelName == ResistanceModel::ONE_D_CONTINUOUS_MODEL) {
        // set correct droplet resistances
        for (auto& [key, channel] : chip->getChannels()) {
            double avgViscosityInChannel = getAverageViscosityInChannel(key);
            double resistance = resistanceModel->getChannelResistance(channel.get(), avgViscosityInChannel);
            channel->setChannelResistance(resistance);
        }
    }
}

void Simulation::updateDropletResistances() {
    // set all droplet resistances of all channels to 0.0
    for (auto& [key, channel] : chip->getChannels()) {
        channel->setDropletResistance(0.0);
    }

    // set correct droplet resistances
    for (auto& [key, droplet] : droplets) {
        // only consider droplets that are inside the network (i.e., also trapped droplets)
        if (droplet->getDropletState() == DropletState::INJECTION || droplet->getDropletState() == DropletState::SINK) {
            continue;
        }

        droplet->addResistances(*resistanceModel);
    }
}

std::vector<std::unique_ptr<Event>> Simulation::computeEvents(double simulationResultTimeCounter) {
    // events
    std::vector<std::unique_ptr<Event>> events;

    // injection events
    for (auto& [key, injection] : injections) {
        double injectionTime = injection->getInjectionTime();
        if (injection->getDroplet()->getDropletState() == DropletState::INJECTION) {
            events.push_back(std::make_unique<InjectionEvent>(injectionTime - currTime, *injection));
        }
    }

    // continuous injection events
    for (auto& [key, changeInputMixture] : changeInputMixtures) {
        if (!changeInputMixture->wasPerformed()) {
            double injectionTime = changeInputMixture->getInjectionTime();
            events.push_back(std::make_unique<ChangeInputMixtureEvent>(injectionTime - currTime, *changeInputMixture));
        }
    }

    if (!changeInputMixtures.empty()) {
        events.push_back(std::make_unique<StopSimulationEvent>(simulationDuration - currTime, *this));
    }

    // define maps that are used for detecting merging inside channels
    std::unordered_map<int, std::vector<DropletBoundary*>> channelBoundariesMap;
    std::unordered_map<DropletBoundary*, Droplet*> boundaryDropletMap;

    for (auto& [key, droplet] : droplets) {
        // only consider droplets inside the network (but no trapped droplets)
        if (droplet->getDropletState() != DropletState::NETWORK) {
            continue;
        }

        // loop through boundaries
        for (auto& boundary : droplet->getBoundaries()) {
            // the flow rate of the boundary indicates if the boundary moves towards or away from the droplet center and, hence, if a BoundaryTailEvent or BoundaryHeadEvent should occur, respectively
            // if the flow rate of the boundary is 0, then no events will be triggered (the boundary may be in a Wait state)
            if (boundary->getFlowRate() < 0) {
                // boundary moves towards the droplet center => BoundaryTailEvent
                double time = boundary->getTime();
                events.push_back(std::make_unique<BoundaryTailEvent>(time, *droplet, *boundary, *chip));
            } else if (boundary->getFlowRate() > 0) {
                // boundary moves away from the droplet center => BoundaryHeadEvent
                double time = boundary->getTime();

                // in this scenario also a MergeBifurcationEvent can happen when merging is enabled
                // this means a boundary comes to a bifurcation where a droplet is already present
                // hence it is either a MergeBifurcationEvent or a BoundaryHeadEvent that will happen

                // check if merging is enabled
                Droplet* mergeDroplet = nullptr;
                if (enableMerging) {
                    // find droplet to merge (if present)
                    auto referenceNode = boundary->getOppositeReferenceNode();
                    mergeDroplet = getDropletAtNode(referenceNode->getId());
                }

                if (mergeDroplet == nullptr) {
                    // no merging will happen => BoundaryHeadEvent
                    events.push_back(std::make_unique<BoundaryHeadEvent>(time, *droplet, *boundary, *chip));
                } else {
                    // merging of the actual droplet with the merge droplet will happen => MergeBifurcationEvent
                    events.push_back(std::make_unique<MergeBifurcationEvent>(time, *droplet, *boundary, *mergeDroplet, *this));
                }
            }

            // fill the maps which are later used for merging inside channels (if merging is enabled)
            if (enableMerging) {
                auto [value, success] = channelBoundariesMap.try_emplace(boundary->getChannelPosition().getChannel()->getId());
                value->second.push_back(boundary.get());

                boundaryDropletMap.emplace(boundary.get(), droplet.get());
            }
        }
    }

    // check for MergeChannelEvents, i.e, for boundaries of other droplets that are in the same channel
    // here the previously defined maps are used => if merging is not enabled these maps are emtpy
    // loop through channelsBoundariesMap
    for (auto& [channelId, boundaries] : channelBoundariesMap) {
        // loop through boundaries that are inside this channel
        for (size_t i = 0; i < boundaries.size(); i++) {
            // get reference boundary and droplet
            auto referenceBoundary = boundaries[i];
            auto referenceDroplet = boundaryDropletMap.at(referenceBoundary);

            // get channel
            auto channel = referenceBoundary->getChannelPosition().getChannel();

            // get velocity and absolute position of the boundary
            // positive values for v0 indicate a movement from node0 towards node1
            auto q0 = referenceBoundary->isVolumeTowardsNode0() ? referenceBoundary->getFlowRate() : -referenceBoundary->getFlowRate();
            auto v0 = q0 / channel->getArea();
            auto p0 = referenceBoundary->getChannelPosition().getPosition() * channel->getLength();

            // compare reference boundary against all others
            // the two loops are defined in such a way, that only one merge event happens for a pair of boundaries
            for (size_t j = i + 1; j < boundaries.size(); j++) {
                auto boundary = boundaries[j];
                auto droplet = boundaryDropletMap.at(boundary);

                // do not consider if this boundary is form the same droplet
                if (droplet == referenceDroplet) {
                    continue;
                }

                // get velocity and absolute position of the boundary
                // positive values for v0 indicate a movement from node0 towards node1
                auto q1 = boundary->isVolumeTowardsNode0() ? boundary->getFlowRate() : -boundary->getFlowRate();
                auto v1 = q1 / channel->getArea();
                auto p1 = boundary->getChannelPosition().getPosition() * channel->getLength();

                // do not merge when both velocities are equal (would result in infinity time)
                if (v0 == v1) {
                    continue;
                }

                // compute time and merge position
                auto time = (p1 - p0) / (v0 - v1);
                auto pMerge = p0 + v0 * time;  // or p1 + v1*time
                auto pMergeRelative = pMerge / channel->getLength();

                // do not trigger a merge event when:
                // * time is negative => indicates that both boundaries go in different directions or that one boundary cannot "outrun" the other because it is too slow
                // * relative merge position is outside the range of [0, 1] => the merging would happen "outside" the channel and a boundary would already switch a channel before this event could happen
                if (time < 0 || pMergeRelative < 0 || 1 < pMergeRelative) {
                    continue;
                }

                // add MergeChannelEvent
                events.push_back(std::make_unique<MergeChannelEvent>(time, *referenceDroplet, *referenceBoundary, *droplet, *boundary, *this));
            }
        }
    }

    // time step event
    if ((dropletsAtBifurcation && maximalAdaptiveTimeStep > 0)) {
        events.push_back(std::make_unique<TimeStepEvent>(maximalAdaptiveTimeStep));
    }

    if (!changeInputMixtures.empty()) {  //only simulate all minimal timesteps if continuous fluid simulation
        events.push_back(std::make_unique<TimeStepEvent>(std::min(simulationResultTimeCounter, internalSimulationTimeStep)));
    }

    return events;
}

void Simulation::moveDroplets(double timeStep) {
    // loop over all droplets
    for (auto& [key, droplet] : droplets) {
        // only consider droplets inside the network (but no trapped droplets)
        if (droplet->getDropletState() != DropletState::NETWORK) {
            continue;
        }

        // loop through boundaries
        for (auto& boundary : droplet->getBoundaries()) {
            // move boundary in correct direction
            boundary->moveBoundary(timeStep);
        }
    }
}

double Simulation::getPumpOutflowVolume(int nodeId, double flowRate, double timeStep) {
    double totalArea = 0.0;
    double totalVolume = 0.0;
    for (auto channel : chip->getChannelsAtNode(nodeId)) {
        totalArea += channel->getArea();
        totalVolume += channel->getVolume();
    }
    double movedDistance = (std::abs(flowRate) * timeStep) / totalVolume;
    double outflowVolume = movedDistance * totalArea;
    return outflowVolume;
}

void Simulation::calculateNewMixtures(double timeStep) {
    /* 
    move positions and calculate mixture inflow at nodes 
    */
    struct MixtureInflow {
        int mixtureId;
        double inflowVolume;
    };
    std::unordered_map<int, std::vector<MixtureInflow>> mixtureInflowAtNode;  // <nodeId <mixtureId, inflowVolume>>
    std::unordered_map<int, double> totalInflowVolumeAtNode;

    for (auto& [nodeId, node] : chip->getNodes()) {
        // pumps
        for (auto& [key, pressurePump] : chip->getPressurePumps()) {
            auto flowRate = pressurePump->getFlowRate();
            auto pumpNodeId = flowRate > 0.0 ? pressurePump->getNode1()->getId() : pressurePump->getNode0()->getId();
            if (pumpNodeId == nodeId) {
                double inflowVolume = getPumpOutflowVolume(nodeId, flowRate, timeStep);
                MixtureInflow mixtureInflow = {pressurePump->getMixtureId(), inflowVolume};
                mixtureInflowAtNode[nodeId].push_back(mixtureInflow);
                auto [iterator, inserted] = totalInflowVolumeAtNode.try_emplace(nodeId, inflowVolume);
                if (!inserted) {
                    iterator->second = iterator->second + inflowVolume;
                }
            }
        }
        for (auto& [key, flowRatePump] : chip->getFlowRatePumps()) {
            auto flowRate = flowRatePump->getFlowRate();
            auto pumpNodeId = flowRate > 0.0 ? flowRatePump->getNode1()->getId() : flowRatePump->getNode0()->getId();
            if (pumpNodeId == nodeId) {
                double inflowVolume = getPumpOutflowVolume(nodeId, flowRate, timeStep);
                MixtureInflow mixtureInflow = {flowRatePump->getMixtureId(), inflowVolume};
                mixtureInflowAtNode[nodeId].push_back(mixtureInflow);
                auto [iterator, inserted] = totalInflowVolumeAtNode.try_emplace(nodeId, inflowVolume);
                if (!inserted) {
                    iterator->second = iterator->second + inflowVolume;
                }
            }
        }

        // channels
        for (auto channel : chip->getChannelsAtNode(nodeId)) {
            // only consider channels where there is an inflow at this node (= where the flow direction is into the node)
            auto flowRate = channel->getFlowRate();
            if ((flowRate > 0.0 && channel->getNode1()->getId() == nodeId) || (flowRate < 0.0 && channel->getNode0()->getId() == nodeId)) {
                for (auto& [mixtureId, endPos] : mixturesInEdge.at(channel->getId())) {
                    double movedDistance = (std::abs(flowRate) * timeStep) / channel->getVolume();
                    double newEndPos = 0.0;
                    double inflowVolume = 0.0;
                    double channelArea = channel->getArea();
                    if (flowRate < 0) {
                        newEndPos = std::max(0.0, endPos - movedDistance);
                        inflowVolume = (newEndPos * -1) * channelArea + movedDistance * channelArea;
                    } else {
                        newEndPos = std::min(endPos + movedDistance, 1.0);
                        inflowVolume = (newEndPos - 1.0) * channelArea + movedDistance * channelArea;
                    }
                    endPos = newEndPos;
                    if ((flowRate > 0 && newEndPos == 1.0) || (flowRate < 0 && newEndPos == 0.0)) {
                        // fluid flows into node, add to mixture inflow
                        MixtureInflow mixtureInflow = {mixtureId, inflowVolume};
                        mixtureInflowAtNode[nodeId].push_back(mixtureInflow);
                        auto [iterator, inserted] = totalInflowVolumeAtNode.try_emplace(nodeId, inflowVolume);
                        if (!inserted) {
                            iterator->second = iterator->second + inflowVolume;
                        }
                    }
                }
                // with a negative flow-rate, the last mixture in the channel spans from pos 0.0 to xy
                // therefore the outlow at the 0 node must also calculated from pos 0.0
                if (flowRate < 0.0 && channel->getNode0()->getId() == nodeId) {
                    auto& [mixtureId, endPos] = mixturesInEdge.at(channel->getId()).back();
                    double movedDistance = (std::abs(flowRate) * timeStep) / channel->getVolume();
                    double newEndPos = 0.0;
                    double inflowVolume = 0.0;
                    double channelArea = channel->getArea();
                    if (flowRate < 0) {
                        newEndPos = movedDistance;
                        inflowVolume = movedDistance * channelArea;
                    }

                    // fluid flows into node, add to mixture inflow
                    MixtureInflow mixtureInflow = {mixtureId, inflowVolume};
                    mixtureInflowAtNode[nodeId].push_back(mixtureInflow);
                    auto [iterator, inserted] = totalInflowVolumeAtNode.try_emplace(nodeId, inflowVolume);
                    if (!inserted) {
                        iterator->second = iterator->second + inflowVolume;
                    }
                }
            }
        }

        // membranes & organs
        // only update fluid positions
        // contributes to outflow at node only indirectly through the fluid concentration exchange with the connected channel
        for (auto& membrane : chip->getMembranesAtNode(nodeId)) {
            auto* organ = membrane->getOrgan();
            auto* channel = membrane->getChannel();
            // mixtures move through the organ at the same speed as through the connected channel
            // this is an abstraction to get time-accurate results
            // in reality, there is no flow rate in the organ
            auto flowRate = channel->getFlowRate();
            if ((flowRate > 0.0 && channel->getNode1()->getId() == nodeId) || (flowRate < 0.0 && channel->getNode0()->getId() == nodeId)) {
                for (auto& [mixtureId, endPos] : mixturesInEdge.at(membrane->getId())) {
                    double movedDistance = (std::abs(flowRate) * timeStep) / channel->getVolume();
                    double newEndPos = 0.0;
                    if (flowRate < 0) {
                        newEndPos = std::max(0.0, endPos - movedDistance);
                    } else {
                        newEndPos = std::min(endPos + movedDistance, 1.0);
                    }
                    endPos = newEndPos;
                }
                for (auto& [mixtureId, endPos] : mixturesInEdge.at(organ->getId())) {
                    double movedDistance = (std::abs(flowRate) * timeStep) / channel->getVolume();
                    double newEndPos = 0.0;
                    if (flowRate < 0) {
                        newEndPos = std::max(0.0, endPos - movedDistance);
                    } else {
                        newEndPos = std::min(endPos + movedDistance, 1.0);
                    }
                    endPos = newEndPos;
                }
            }
        }
    }

    /*
     calculate mixture outflow at node from inflow
     */
    std::unordered_map<int, int> mixtureOutflowAtNode;
    for (auto& [nodeId, mixtureInflowList] : mixtureInflowAtNode) {
        std::unordered_map<int, double> newFluidConcentrations;
        for (auto& mixtureInflow : mixtureInflowList) {
            for (auto& [fluidId, oldConcentration] : mixtures.at(mixtureInflow.mixtureId).getFluidConcentrations()) {
                double newConcentration = oldConcentration * mixtureInflow.inflowVolume / totalInflowVolumeAtNode.at(nodeId);
                auto [iterator, inserted] = newFluidConcentrations.try_emplace(fluidId, newConcentration);
                if (!inserted) {
                    iterator->second = iterator->second + newConcentration;
                }
            }
        }

        int newMixtureId = addMixture(std::move(newFluidConcentrations));
        mixtureOutflowAtNode.try_emplace(nodeId, newMixtureId);
    }

    /* 
    add outflow as inflow to edges
    */
    for (auto& [nodeId, node] : chip->getNodes()) {
        // channels
        for (auto& channel : chip->getChannelsAtNode(nodeId)) {
            // check if edge is an outflow edge to this node
            auto flowRate = channel->getFlowRate();
            if ((flowRate > 0.0 && channel->getNode0()->getId() == nodeId) || (flowRate < 0.0 && channel->getNode1()->getId() == nodeId)) {
                double newPos = std::abs(flowRate) * timeStep / channel->getVolume();
                assert(newPos <= 1.0 && newPos >= 0.0);

                bool oldEqualsNewConcentration = true;
                auto& oldFluidConcentrations = mixtures.at(mixturesInEdge.at(channel->getId()).back().first).getFluidConcentrations();
                if (mixtureOutflowAtNode.count(nodeId)) {
                    if (flowRate < 0.0) {
                        mixturesInEdge.at(channel->getId()).push_front(std::make_pair(mixtureOutflowAtNode.at(nodeId), 1.0));
                    } else {
                        mixturesInEdge.at(channel->getId()).push_back(std::make_pair(mixtureOutflowAtNode.at(nodeId), newPos));
                    }
                }
            }
        }

        // membranes & organs
        for (auto& membrane : chip->getMembranesAtNode(nodeId)) {
            auto* organ = membrane->getOrgan();
            auto* membraneChannel = membrane->getChannel();
            double channelFlowRate = membrane->getChannel()->getFlowRate();
            // check if edge is an outflow edge to this node
            if ((channelFlowRate > 0.0 && organ->getNode0()->getId() == nodeId) || (channelFlowRate < 0.0 && organ->getNode1()->getId() == nodeId)) {
                double movedDistance = (std::abs(channelFlowRate) * timeStep) / membraneChannel->getVolume();
                double newEndPos = movedDistance;
                assert(newEndPos <= 1.0 && newEndPos >= 0.0);
                // copy the mixture that outflows the organ as new mixture at the inflow (to ensure mass conservation)
                // concentration that diffuses through membrane from new inflow is added later (see below)
                auto& currMixtureOrgan = mixtures.at(mixturesInEdge.at(organ->getId()).front().first);
                std::unordered_map<int, double> newFluidConcentrations(currMixtureOrgan.getFluidConcentrations());
                int newMixtureId = addMixture(std::move(newFluidConcentrations));

                if (mixtureOutflowAtNode.count(nodeId)) {
                    if (channelFlowRate < 0.0) {
                        mixturesInEdge.at(membrane->getId()).push_front(std::make_pair(newMixtureId, 1.0));
                        mixturesInEdge.at(organ->getId()).push_front(std::make_pair(newMixtureId, 1.0));
                    } else {
                        mixturesInEdge.at(membrane->getId()).push_back(std::make_pair(newMixtureId, newEndPos));
                        mixturesInEdge.at(organ->getId()).push_back(std::make_pair(newMixtureId, newEndPos));
                    }
                }
            }
        }
    }

    /*
    calculate exchange between organ and channel through membranes and change mixtures accordingly
    */
    for (auto& [nodeId, node] : chip->getNodes()) {
        for (auto membrane : chip->getMembranesAtNode(nodeId)) {
            auto* organ = membrane->getOrgan();
            // mixtures move through the organ at the same speed as through the connected channel
            // this is an abstraction to get time-accurate results
            // in reality, there is no flow rate in the organ
            auto* channel = membrane->getChannel();
            double channelFlowRate = channel->getFlowRate();
            if ((channelFlowRate > 0.0 && organ->getNode1()->getId() == nodeId) || (channelFlowRate < 0.0 && organ->getNode0()->getId() == nodeId)) {
                int mixturesInChannelSize = mixturesInEdge.at(channel->getId()).size();
                int mixturesInOrganSize = mixturesInEdge.at(organ->getId()).size();
                assert(mixturesInChannelSize == mixturesInOrganSize);
                int dequeIdx = mixturesInChannelSize - 1;
                double startPos = 0.0;
                for (auto it = mixturesInEdge.at(organ->getId()).rbegin(); it != mixturesInEdge.at(organ->getId()).rend(); it++) {
                    auto& [oldOrganMixtureId, endPos] = *it;
                    double mixtureLengthAbs = (endPos - startPos) * channel->getLength();

                    auto& currMixtureOrgan = mixtures.at(oldOrganMixtureId);
                    std::unordered_map<int, double> newFluidConcentrationsOrgan(currMixtureOrgan.getFluidConcentrations());
                    int newOrganMixtureId = addMixture(std::move(newFluidConcentrationsOrgan));
                    oldOrganMixtureId = newOrganMixtureId;

                    auto& currMixtureChannel = mixtures.at(mixturesInEdge.at(channel->getId()).at(dequeIdx).first);
                    std::unordered_map<int, double> newFluidConcentrationsChannel(currMixtureChannel.getFluidConcentrations());
                    int newChannelMixtureId = addMixture(std::move(newFluidConcentrationsChannel));
                    mixturesInEdge.at(channel->getId()).at(dequeIdx).first = newChannelMixtureId;

                    for (auto& [fluidId, fluid] : fluids) {
                        double area = membrane->getWidth() * mixtureLengthAbs;
                        double resistance = membraneResistanceModel->getMembraneResistance(membrane, fluids.at(fluidId).get(), area);
                        double fluidSaturation = fluids.at(fluidId)->getSaturation();
                        if (fluidSaturation != 0.0 && mixtureLengthAbs > 0.0) {
                            double concentrationChannel = mixtures.at(newChannelMixtureId).getConcentrationOfFluid(fluidId);
                            double concentrationOrgan = mixtures.at(newOrganMixtureId).getConcentrationOfFluid(fluidId);

                            // positive flux defined to go from channel to organ
                            // negative flux defined to go from organ to channel
                            double concentrationDifference = (concentrationChannel - concentrationOrgan);
                            double concentrationChangeMol = membrane->getConcentrationChange(resistance, timeStep, concentrationDifference, currTime);

                            double concentrationChangeOrgan = concentrationChangeMol / (mixtureLengthAbs * organ->getWidth() * organ->getHeight());
                            double concentrationChangeChannel = concentrationChangeMol * -1 / (mixtureLengthAbs * channel->getWidth() * channel->getHeight());
                            mixtures.at(newOrganMixtureId).changeFluidConcentration(fluidId, concentrationChangeOrgan);
                            mixtures.at(newChannelMixtureId).changeFluidConcentration(fluidId, concentrationChangeChannel);
                        }
                    }
                    startPos = endPos;
                    dequeIdx--;
                }
            }
        }
    }

    /*
    remove fluids that have outflowed their edge
    */
    for (auto& [nodeId, node] : chip->getNodes()) {
        // channels
        for (auto& channel : chip->getChannelsAtNode(nodeId)) {
            auto count = 0;  // to not remove the 1.0 fluid if there is only one
            for (auto& [mixtureId, endPos] : mixturesInEdge.at(channel->getId())) {
                //remove mixtures that completely flow out of channel (only 1 fluid with pos 1.0 or 0.0 left)
                if (endPos == 0.0) {
                    if (count != 0) {
                        mixturesInEdge.at(channel->getId()).pop_back();
                    }
                    count++;
                }
                if (endPos == 1.0) {
                    if (count != 0) {
                        mixturesInEdge.at(channel->getId()).pop_front();
                    }
                    count++;
                }
            }
        }
        // membranes and organs
        for (auto& membrane : chip->getMembranesAtNode(nodeId)) {
            auto count = 0;  // to not remove the 1.0 fluid if there is only one
            auto* organ = membrane->getOrgan();
            auto* channel = membrane->getChannel();

            for (auto& [mixtureId, endPos] : mixturesInEdge.at(membrane->getId())) {
                //count mixtures that completely flow out of channel (only 1 fluid with pos 1.0 left)
                if (endPos == 0.0) {
                    if (count != 0) {
                        mixturesInEdge.at(membrane->getId()).pop_back();
                        mixturesInEdge.at(membrane->getOrgan()->getId()).pop_back();
                    }
                    count++;
                }
                if (endPos == 1.0) {
                    if (count != 0) {
                        mixturesInEdge.at(membrane->getId()).pop_front();
                        mixturesInEdge.at(membrane->getOrgan()->getId()).pop_front();
                    }
                    count++;
                }
            }
        }
    }
}

void Simulation::storeSimulationParameters(droplet::SimulationResult& result) {
    // chip name
    result.chip.name = chip->getName();

    // channel
    for (const auto& [key, channel] : chip->getChannels()) {
        result.chip.channels.try_emplace(key, channel->getId(), channel->getName(), channel->getNode0()->getId(), channel->getNode1()->getId(), channel->getWidth(), channel->getHeight(), channel->getLength(), static_cast<droplet::ChannelType>(static_cast<int>(channel->getChannelType())));
    }

    // membrane
    for (const auto& [key, membrane] : chip->getMembranes()) {
        result.chip.membranes.try_emplace(key, membrane->getId(), membrane->getName(), membrane->getNode0()->getId(), membrane->getNode1()->getId(), membrane->getWidth(), membrane->getHeight(), membrane->getLength(), membrane->getPoreRadius(), membrane->getPorosity(), membrane->getChannel()->getId(), membrane->getOrgan()->getId());
    }

    // organ
    for (const auto& [key, organ] : chip->getOrgans()) {
        result.chip.organs.try_emplace(key, organ->getId(), organ->getName(), organ->getNode0()->getId(), organ->getNode1()->getId(), organ->getWidth(), organ->getHeight(), organ->getLength());
    }

    // flow rate pump
    for (const auto& [key, pump] : chip->getFlowRatePumps()) {
        result.chip.flowRatePumps.try_emplace(key, pump->getId(), pump->getName(), pump->getNode0()->getId(), pump->getNode1()->getId(), pump->getFlowRate());
    }

    // pressure pump
    for (const auto& [key, pump] : chip->getPressurePumps()) {
        result.chip.pressurePumps.try_emplace(key, pump->getId(), pump->getName(), pump->getNode0()->getId(), pump->getNode1()->getId(), pump->getPressure());
    }

    // fluids
    for (const auto& [key, fluid] : fluids) {
        auto [value, success] = result.fluids.try_emplace(key, fluid->getId(), fluid->getName(), fluid->getViscosity(), fluid->getDensity(), fluid->getConcentration(), fluid->getMolecularSize(), fluid->getDiffusionCoefficient(), fluid->getSaturation());

        //mixed fluids
        for (const auto mixedFluid : fluid->getMixedFluids()) {
            value->second.mixedFluidIds.push_back(mixedFluid->getId());
        }
    }

    // droplet
    for (const auto& [key, droplet] : droplets) {
        auto [value, success] = result.droplets.try_emplace(key, droplet->getId(), droplet->getName(), droplet->getVolume(), droplet->getFluid()->getId());

        //merged droplets
        for (const auto mergedDroplet : droplet->getMergedDroplets()) {
            value->second.mergedDropletIds.push_back(mergedDroplet->getId());
        }
    }

    // injections
    for (const auto& [key, injection] : injections) {
        result.injections.try_emplace(key, injection->getId(), injection->getDroplet()->getId(), injection->getInjectionTime(), injection->getInjectionPosition().getChannel()->getId(), injection->getInjectionPosition().getPosition());
    }

    for (const auto& [key, continuousInjection] : changeInputMixtures) {
        result.continuousInjections.try_emplace(key, continuousInjection->getId(), continuousInjection->getMixtureId(), continuousInjection->getInjectionTime(), continuousInjection->getInjectionPump()->getId());
    }

    result.continuousPhaseMixtureId = continuousPhaseMixtureId;
    result.maximalAdaptiveTimeStep = maximalAdaptiveTimeStep;

    switch (resistanceModelName) {
        case ResistanceModel::ONE_D_MODEL:
            result.resistanceModel = 0;
            break;
        case ResistanceModel::TEST_MODEL:
            result.resistanceModel = 1;
            break;
        case ResistanceModel::ONE_D_CONTINUOUS_MODEL:
            result.resistanceModel = 2;
            break;
    }

    switch (membraneResistanceModelName) {
        case MembraneResistanceModel::MODEL_0:
            result.membraneResistanceModel = 0;
            break;
        case MembraneResistanceModel::MODEL_1:
            result.membraneResistanceModel = 1;
            break;
        case MembraneResistanceModel::MODEL_2:
            result.membraneResistanceModel = 2;
            break;
        case MembraneResistanceModel::MODEL_3:
            result.membraneResistanceModel = 3;
            break;
        case MembraneResistanceModel::MODEL_4:
            result.membraneResistanceModel = 4;
            break;
        case MembraneResistanceModel::MODEL_5:
            result.membraneResistanceModel = 5;
            break;
        case MembraneResistanceModel::MODEL_6:
            result.membraneResistanceModel = 6;
            break;
        case MembraneResistanceModel::MODEL_7:
            result.membraneResistanceModel = 7;
            break;
        case MembraneResistanceModel::MODEL_8:
            result.membraneResistanceModel = 8;
            break;
        case MembraneResistanceModel::MODEL_9:
            result.membraneResistanceModel = 9;
            break;
    }
}

void Simulation::storeSimulationResults(droplet::SimulationResult& result) {
    // add new fluids if present (try_emplace does nothing if key is already present)
    for (const auto& [key, fluid] : fluids) {
        auto [value, success] = result.fluids.try_emplace(key, fluid->getId(), fluid->getName(), fluid->getViscosity(), fluid->getDensity(), fluid->getConcentration(), fluid->getMolecularSize(), fluid->getDiffusionCoefficient(), fluid->getSaturation());

        //mixed fluids
        if (success) {
            for (const auto mixedFluid : fluid->getMixedFluids()) {
                value->second.mixedFluidIds.push_back(mixedFluid->getId());
            }
        }

        //store inital mixtures that consist of 100% of one fluid
        auto& mixture = mixtures.at(fluid->getMixtureId());
        result.mixtures.try_emplace(mixture.getId(), mixture.getId(), mixture.getFluidConcentrations(), mixture.getViscosity(), mixture.getDensity(), mixture.getLargestMolecularSize());
    }

    // add new droplets if present (try_emplace does nothing if key is already present)
    for (const auto& [key, droplet] : droplets) {
        auto [value, success] = result.droplets.try_emplace(key, droplet->getId(), droplet->getName(), droplet->getVolume(), droplet->getFluid()->getId());

        //merged droplets
        if (success) {
            for (const auto mergedDroplet : droplet->getMergedDroplets()) {
                value->second.mergedDropletIds.push_back(mergedDroplet->getId());
            }
        }
    }

    // state
    auto& state = result.states.emplace_back(iState++, currTime);

    // pressures
    for (auto& [id, node] : chip->getNodes()) {
        state.pressures.try_emplace(id, node->getPressure());
    }

    // flow rates
    for (auto& [id, channel] : chip->getChannels()) {
        state.flowRates.try_emplace(id, channel->getFlowRate());
    }
    for (auto& [id, pump] : chip->getFlowRatePumps()) {
        state.flowRates.try_emplace(id, pump->getFlowRate());
    }
    for (auto& [id, pump] : chip->getPressurePumps()) {
        state.flowRates.try_emplace(id, pump->getFlowRate());
    }

    // droplet positions
    for (auto& [id, droplet] : droplets) {
        // create new droplet position
        auto [value, success] = state.dropletPositions.try_emplace(id, static_cast<droplet::DropletState>(static_cast<int>(droplet->getDropletState())));

        // add boundaries
        for (auto& boundary : droplet->getBoundaries()) {
            // get channel position
            auto channelPosition = boundary->getChannelPosition();

            // add boundary
            value->second.boundaries.emplace_back(channelPosition.getChannel()->getId(), channelPosition.getPosition(), boundary->isVolumeTowardsNode0(), boundary->getFlowRate(), static_cast<droplet::BoundaryState>(static_cast<int>(boundary->getState())));
        }

        // add fully occupied channels
        for (auto& channel : droplet->getFullyOccupiedChannels()) {
            // add channel
            value->second.channelIds.emplace_back(channel->getId());
        }
    }

    //fluidInEdge
    for (auto& [edgeId, elem] : mixturesInEdge) {
        std::deque<std::pair<int, double>> mixturePos = elem;
        state.mixturesInEdge.try_emplace(edgeId, mixturePos);

        //only store mixtures that are used
        for (auto& [mixtureId, _mixturePos] : mixturePos) {
            auto& mixture = mixtures.at(mixtureId);
            result.mixtures.try_emplace(mixture.getId(), mixture.getId(), mixture.getFluidConcentrations(), mixture.getViscosity(), mixture.getDensity(), mixture.getLargestMolecularSize());
        }
    }
}

Fluid* Simulation::mixFluids(int fluid0Id, double volume0, int fluid1Id, double volume1) {
    // check if fluids are identically (no merging needed) and if they exist
    if (fluid0Id == fluid1Id) {
        // try to get the fluid (throws error if the fluid is not present)
        return fluids.at(fluid0Id).get();
    }

    // get fluids
    auto fluid0 = fluids.at(fluid0Id).get();
    auto fluid1 = fluids.at(fluid1Id).get();

    // compute ratios
    double volume = volume0 + volume1;
    double ratio0 = volume0 / volume;
    double ratio1 = volume1 / volume;

    // compute new fluid values
    double viscosity = ratio0 * fluid0->getViscosity() + ratio1 * fluid1->getViscosity();
    double density = ratio0 * fluid0->getDensity() + ratio1 * fluid1->getDensity();
    double concentration = ratio0 * fluid0->getConcentration() + ratio1 * fluid1->getConcentration();

    // add new fluid
    auto newFluid = addFluid(viscosity, density, concentration, fluid0->getMolecularSize(), fluid0->getDiffusionCoefficient(), fluid0->getSaturation());

    //add previous fluids
    newFluid->addMixedFluid(fluid0);
    newFluid->addMixedFluid(fluid1);

    return newFluid;
}

Droplet* Simulation::mergeDroplets(int droplet0Id, int droplet1Id) {
    // check if droplets are identically (no merging needed) and if they exist
    if (droplet0Id == droplet1Id) {
        // try to get the droplet (throws error if the droplet is not present)
        return droplets.at(droplet0Id).get();
    }

    // get droplets
    auto droplet0 = getDroplet(droplet0Id);
    auto droplet1 = getDroplet(droplet1Id);

    // compute volumes
    double volume0 = droplet0->getVolume();
    double volume1 = droplet1->getVolume();
    double volume = volume0 + volume1;

    // merge fluids
    auto fluid = mixFluids(droplet0->getFluid()->getId(), volume0, droplet1->getFluid()->getId(), volume1);

    // add new droplet
    auto newDroplet = addDroplet(fluid->getId(), volume);

    //add previous droplets
    newDroplet->addMergedDroplet(droplet0);
    newDroplet->addMergedDroplet(droplet1);

    return newDroplet;
}

void Simulation::stopSimulation() {
    simulationInProgress = false;
}

}  // namespace sim