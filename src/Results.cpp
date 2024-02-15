#include "Results.h"

#include <deque>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "nlohmann/json.hpp"

namespace droplet {

ChannelPosition::ChannelPosition(int channelId, double position) : channelId(channelId), position(position) {}

DropletBoundary::DropletBoundary(int channelId, double position, bool volumeTowards0, double flowRate, BoundaryState state) : position(ChannelPosition(channelId, position)), volumeTowards0(volumeTowards0), flowRate(flowRate), state(state) {}

DropletPosition::DropletPosition(DropletState state) : state(state) {}

Fluid::Fluid(int id, std::string name, double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient, double saturation) : id(id), name(name), viscosity(viscosity), density(density), concentration(concentration), molecularSize(molecularSize), diffusionCoefficient(diffusionCoefficient), saturation(saturation) {}

Mixture::Mixture(int id, std::unordered_map<int, double> fluidConcentrations, double viscosity, double density, double largestMolecularSize) : id(id), fluidConcentrations(fluidConcentrations), viscosity(viscosity), density(density), largestMolecularSize(largestMolecularSize) {}

Droplet::Droplet(int id, std::string name, double volume, int fluidId) : id(id), name(name), volume(volume), fluidId(fluidId) {}

Injection::Injection(int id, int dropletId, double time, int channelId, double position) : id(id), dropletId(dropletId), time(time), position(channelId, position) {}

ContinuousInjection::ContinuousInjection(int id, int mixtureId, double time, int pumpId) : id(id), mixtureId(mixtureId), time(time), pumpId(pumpId) {}

Channel::Channel(int id, std::string name, int node0Id, int node1Id, double width, double height, double length, ChannelType type) : id(id), name(name), node0Id(node0Id), node1Id(node1Id), width(width), height(height), length(length), type(type) {}

Membrane::Membrane(int id, std::string name, int node0Id, int node1Id, double width, double height, double length, double poreRadius, double porosity, int channelId, int organId) : id(id), name(name), node0Id(node0Id), node1Id(node1Id), width(width), height(height), length(length), poreRadius(poreRadius), porosity(porosity), channelId(channelId), organId(organId) {}

Organ::Organ(int id, std::string name, int node0Id, int node1Id, double width, double height, double length) : id(id), name(name), node0Id(node0Id), node1Id(node1Id), width(width), height(height), length(length) {}

FlowRatePump::FlowRatePump(int id, std::string name, int node0Id, int node1Id, double flowRate) : id(id), name(name), node0Id(node0Id), node1Id(node1Id), flowRate(flowRate) {}

PressurePump::PressurePump(int id, std::string name, int node0Id, int node1Id, double pressure) : id(id), name(name), node0Id(node0Id), node1Id(node1Id), pressure(pressure) {}

State::State(int id, double time) : id(id), time(time) {}

double State::getPressure(int nodeId) const {
    return pressures.at(nodeId);
}

double State::getPressureDrop(int node0Id, int node1Id) const {
    return getPressure(node0Id) - getPressure(node1Id);
}

double State::getFlowRate(int channelId) const {
    return flowRates.at(channelId);
}

std::deque<std::pair<int, double>> State::getMixturesInEdge(int edgeId) const {
    return mixturesInEdge.at(edgeId);
}

DropletPathPosition::DropletPathPosition(int stateId) : stateId(stateId) {}

DropletPath::DropletPath(int dropletId) : dropletId(dropletId) {}

std::string DropletPath::toJson(int indent) const {
    nlohmann::json json;

    json["dropletId"] = dropletId;
    json["positions"] = nlohmann::json::array();
    for (auto& position : positions) {
        //channelIds
        auto channelIds = nlohmann::json::array();
        for (auto channelId : position.channelIds) {
            channelIds.push_back(channelId);
        }

        json["positions"].push_back({{"stateId", position.stateId}, {"channelIds", channelIds}});
    }

    return json.dump(indent);
}

DropletPath SimulationResult::getDropletPath(int dropletId) const {
    DropletPath dropletPath(dropletId);

    //loop through states
    for (auto const& state : states) {
        //get dropletPosition
        auto itDropletPosition = state.dropletPositions.find(dropletId);
        if (itDropletPosition != state.dropletPositions.end()) {
            //only consider droplets that are inside the Network
            if (itDropletPosition->second.state != DropletState::NETWORK) {
                continue;
            }

            //create DropletPathPosition
            auto& position = dropletPath.positions.emplace_back(state.id);

            //add channelIds of boundaries
            for (auto& boundary : itDropletPosition->second.boundaries) {
                position.channelIds.insert(boundary.position.channelId);
            }

            //add fully occupied channelIds
            for (auto& channelId : itDropletPosition->second.channelIds) {
                position.channelIds.insert(channelId);
            }

            //check if the set with the channelIds is the same as in the previous state
            //if yes then do not consider this DropletPathPosition and pop it from the list (prevents duplicates)
            if (dropletPath.positions.size() > 1 && dropletPath.positions[dropletPath.positions.size() - 2].channelIds == position.channelIds) {
                dropletPath.positions.pop_back();
            }
        }
    }

    return dropletPath;
}

std::unordered_map<int, double> SimulationResult::getAverageFluidConcentrationsInEdge(int stateId, int edgeId) const {
    std::unordered_map<int, double> fluidConcentrations;
    double prevMixturePos = 0.0;
    for (auto it = states.at(stateId).mixturesInEdge.at(edgeId).crbegin(); it != states.at(stateId).mixturesInEdge.at(edgeId).crend(); it++) {
        auto& [mixtureId, mixturePos] = *it;
        for (auto& [fluidId, concentration] : mixtures.at(mixtureId).fluidConcentrations) {
            double newConcentration = concentration * (mixturePos - prevMixturePos);
            auto [iterator, inserted] = fluidConcentrations.try_emplace(fluidId, newConcentration);
            if (!inserted) {
                iterator->second = iterator->second + newConcentration;
            }
        }
        prevMixturePos = mixturePos;
    }
    return fluidConcentrations;
}

std::string SimulationResult::toJson(int indent) const {
    nlohmann::json json;

    json["continuousPhaseMixtureId"] = continuousPhaseMixtureId;
    json["maximalAdaptiveTimeStep"] = maximalAdaptiveTimeStep;
    json["resistanceModel"] = resistanceModel;
    json["membraneResistanceModel"] = membraneResistanceModel;

    // chip
    json["chip"] = nlohmann::json::object();
    json["chip"]["name"] = chip.name;
    json["chip"]["channels"] = nlohmann::json::array();
    for (auto const& [key, channel] : chip.channels) {
        json["chip"]["channels"].push_back({{"id", channel.id}, {"name", channel.name}, {"node0Id", channel.node0Id}, {"node1Id", channel.node1Id}, {"width", channel.width}, {"height", channel.height}, {"length", channel.length}, {"type", channel.type}});
    }
    json["chip"]["membranes"] = nlohmann::json::array();
    for (auto const& [key, membrane] : chip.membranes) {
        json["chip"]["membranes"].push_back({{"id", membrane.id}, {"name", membrane.name}, {"node0Id", membrane.node0Id}, {"node1Id", membrane.node1Id}, {"width", membrane.width}, {"height", membrane.height}, {"length", membrane.length}, {"poreRadius", membrane.poreRadius}, {"porosity", membrane.porosity}, {"channelId", membrane.channelId}, {"organId", membrane.organId}});
    }
    json["chip"]["organs"] = nlohmann::json::array();
    for (auto const& [key, organ] : chip.organs) {
        json["chip"]["organs"].push_back({{"id", organ.id}, {"name", organ.name}, {"node0Id", organ.node0Id}, {"node1Id", organ.node1Id}, {"width", organ.width}, {"height", organ.height}, {"length", organ.length}});
    }
    json["chip"]["flowRatePumps"] = nlohmann::json::array();
    for (auto const& [key, flowRatePump] : chip.flowRatePumps) {
        json["chip"]["flowRatePumps"].push_back({{"id", flowRatePump.id}, {"name", flowRatePump.name}, {"node0Id", flowRatePump.node0Id}, {"node1Id", flowRatePump.node1Id}, {"flowRate", flowRatePump.flowRate}});
    }
    json["chip"]["pressurePumps"] = nlohmann::json::array();
    for (auto const& [key, pressurePump] : chip.pressurePumps) {
        json["chip"]["pressurePumps"].push_back({{"id", pressurePump.id}, {"name", pressurePump.name}, {"node0Id", pressurePump.node0Id}, {"node1Id", pressurePump.node1Id}, {"pressure", pressurePump.pressure}});
    }

    // fluids
    json["fluids"] = nlohmann::json::array();
    for (auto const& [key, fluid] : fluids) {
        json["fluids"].push_back({{"id", fluid.id}, {"name", fluid.name}, {"mixedFluidIds", fluid.mixedFluidIds}, {"viscosity", fluid.viscosity}, {"density", fluid.density}, {"concentration", fluid.concentration}, {"molecularSize", fluid.molecularSize}, {"diffusionCoefficient", fluid.diffusionCoefficient}, {"saturation", fluid.saturation}});
    }

    // mixtures
    json["mixtures"] = nlohmann::json::array();
    for (auto const& [key, mixture] : mixtures) {
        json["mixtures"].push_back({{"id", mixture.id}, {"fluidConcentrations", mixture.fluidConcentrations}, {"viscosity", mixture.viscosity}, {"density", mixture.density}, {"largestMolecularSize", mixture.largestMolecularSize}});
    }

    // droplets
    json["droplets"] = nlohmann::json::array();
    for (auto const& [key, droplet] : droplets) {
        json["droplets"].push_back({{"id", droplet.id}, {"name", droplet.name}, {"mergedDropletIds", droplet.mergedDropletIds}, {"volume", droplet.volume}, {"fluidId", droplet.fluidId}});
    }

    // injections
    json["injections"] = nlohmann::json::array();
    for (auto const& [key, injection] : injections) {
        json["injections"].push_back({{"id", injection.id}, {"dropletId", injection.dropletId}, {"time", injection.time}, {"position", {{"channelId", injection.position.channelId}, {"position", injection.position.position}}}});
    }

    json["continuousInjections"] = nlohmann::json::array();
    for (auto const& [key, continuousInjection] : continuousInjections) {
        json["continuousInjections"].push_back({{"id", continuousInjection.id}, {"mixtureId", continuousInjection.mixtureId}, {"time", continuousInjection.time}, {"pumpId", continuousInjection.pumpId}});
    }

    // states
    json["states"] = nlohmann::json::array();
    for (auto const& state : states) {
        //pressures
        auto pressures = nlohmann::json::object();
        for (auto const& [key, pressure] : state.pressures) {
            pressures[std::to_string(key)] = pressure;
        }

        //flowRates
        auto flowRates = nlohmann::json::object();
        for (auto const& [key, flowRate] : state.flowRates) {
            flowRates[std::to_string(key)] = flowRate;
        }

        //dropletPositions
        auto dropletPositions = nlohmann::json::object();
        for (auto const& [key, dropletPosition] : state.dropletPositions) {
            //dropletPosition
            auto jsonDropletPosition = nlohmann::json::object();

            //state
            jsonDropletPosition["state"] = dropletPosition.state;

            //boundaries
            // clang-format off
            jsonDropletPosition["boundaries"] = nlohmann::json::array();
            for(auto const &boundary : dropletPosition.boundaries) {
                jsonDropletPosition["boundaries"].push_back({
                    {"volumeTowards0", boundary.volumeTowards0},
                    {"flowRate", boundary.flowRate},
                    {"state", boundary.state},
                    {"position", {
                        {"channelId", boundary.position.channelId},
                        {"position", boundary.position.position}}
                    }
                });
            }
            // clang-format on

            //channelIds
            jsonDropletPosition["channelIds"] = nlohmann::json::array();
            for (auto const& channelId : dropletPosition.channelIds) {
                jsonDropletPosition["channelIds"].push_back(channelId);
            }

            //add dropletPosition to array
            dropletPositions[std::to_string(key)] = jsonDropletPosition;
        }

        // mixtureInEdge
        auto mixturesInEdge = nlohmann::json::object();
        for (auto const& [channelId, fluidIds] : state.mixturesInEdge) {
            mixturesInEdge[std::to_string(channelId)] = fluidIds;
        }

        //add state to array
        json["states"].push_back({{"id", state.id}, {"time", state.time}, {"pressures", pressures}, {"flowRates", flowRates}, {"dropletPositions", dropletPositions}, {"mixtureInEdge", mixturesInEdge}});
    }

    return json.dump(indent);
}

SimulationResult SimulationResult::fromJson(std::string jsonString) {
    auto json = nlohmann::json::parse(jsonString);

    SimulationResult results;

    results.continuousPhaseMixtureId = json["continuousPhaseMixtureId"];
    results.maximalAdaptiveTimeStep = json["maximalAdaptiveTimeStep"];
    results.resistanceModel = json["resistanceModel"];

    //###chip###
    //name
    results.chip.name = json["chip"]["name"];

    //channels
    for (auto& channel : json["chip"]["channels"]) {
        results.chip.channels.try_emplace(channel["id"], channel["id"], channel["name"], channel["node0Id"], channel["node1Id"], channel["width"], channel["height"], channel["length"], static_cast<ChannelType>(channel["type"]));
    }

    //flowRatePumps
    for (auto& pump : json["chip"]["flowRatePumps"]) {
        results.chip.flowRatePumps.try_emplace(pump["id"], pump["id"], pump["name"], pump["node0Id"], pump["node1Id"], pump["flowRate"]);
    }

    //pressurePumps
    for (auto& pump : json["chip"]["pressurePumps"]) {
        results.chip.pressurePumps.try_emplace(pump["id"], pump["id"], pump["name"], pump["node0Id"], pump["node1Id"], pump["pressure"]);
    }

    //##fluids###
    for (auto& fluid : json["fluids"]) {
        auto [value, success] = results.fluids.try_emplace(fluid["id"], fluid["id"], fluid["name"], fluid["viscosity"], fluid["density"], fluid["molecularSize"], fluid["concentration"], fluid["diffusionCoefficient"], fluid["saturation"]);
        for (auto& fluidId : fluid["mixedFluidIds"]) {
            value->second.mixedFluidIds.emplace_back(fluidId);
        }
    }

    //##droplets###
    for (auto& droplet : json["droplets"]) {
        auto [value, success] = results.droplets.try_emplace(droplet["id"], droplet["id"], droplet["name"], droplet["volume"], droplet["fluidId"]);
        for (auto& dropletId : droplet["mergedDropletIds"]) {
            value->second.mergedDropletIds.emplace_back(dropletId);
        }
    }

    //##injections###
    for (auto& injection : json["injections"]) {
        results.injections.try_emplace(injection["id"], injection["id"], injection["dropletId"], injection["time"], injection["position"]["channelId"], injection["position"]["position"]);
    }
    for (auto& continuousInjection : json["continuousInjections"]) {
        results.continuousInjections.try_emplace(continuousInjection["id"], continuousInjection["id"], continuousInjection["mixtureId"], continuousInjection["time"], continuousInjection["pumpId"]);
    }

    //###states###
    for (auto& jsonState : json["states"]) {
        //create state (with corresponding id and time)
        auto& state = results.states.emplace_back(jsonState["id"], jsonState["time"]);

        //pressures
        for (auto& [key, pressure] : jsonState["pressures"].items()) {
            state.pressures.try_emplace(std::stoi(key), pressure);
        }

        //flowRates
        for (auto& [key, flowRate] : jsonState["flowRates"].items()) {
            state.flowRates.try_emplace(std::stoi(key), flowRate);
        }

        //dropletPositions
        for (auto& [key, jsonDropletPosition] : jsonState["dropletPositions"].items()) {
            //create dropletPosition
            auto [value, success] = state.dropletPositions.try_emplace(std::stoi(key), static_cast<DropletPosition>(jsonDropletPosition["state"]));

            //boundaries
            for (auto& boundary : jsonDropletPosition["boundaries"]) {
                value->second.boundaries.emplace_back(boundary["position"]["channelId"], boundary["position"]["position"], boundary["volumeTowards0"], boundary["flowRate"], static_cast<BoundaryState>(jsonDropletPosition["state"]));
            }

            //channelIds
            for (auto& channelId : jsonDropletPosition["channelIds"]) {
                value->second.channelIds.emplace_back(channelId);
            }
        }

        //mixtureInEgde
        for (auto& [channelId, mixtureId] : jsonState["mixturesInEdge"].items()) {
            state.mixturesInEdge.try_emplace(std::stoi(channelId), mixtureId);
        }
    }

    return results;
}

}  // namespace droplet