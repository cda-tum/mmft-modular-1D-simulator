#include "Chip.h"

#include <memory>
#include <unordered_map>
#include <stdexcept>

#include "FlowRatePump.h"
#include "PressurePump.h"

namespace arch {

Chip::Chip() {}

Node* Chip::getOrAddNode(int nodeId) {
    auto result = nodes.insert({nodeId, std::make_unique<Node>(nodeId)});

    if (result.second) {
        // insertion happened and we have to add an additional entry into the network
        network.insert_or_assign(nodeId, std::vector<Channel*>{});
    }

    // return raw pointer to the node
    return result.first->second.get();
}

void Chip::visitNodes(int id, std::unordered_map<int, bool>& visitedNodes, std::unordered_map<int, bool>& visitedChannels) {
    const auto net = network.at(id);
    for (Channel* channel : net) {
        if (!(channel->getChannelType() == ChannelType::CLOGGABLE)) {
            visitedNodes.at(id) = true;
            if (visitedChannels[channel->getId()] == false) {
                visitedChannels.at(channel->getId()) = true;
                if (channel->getNode0()->getId() != id) {
                    visitNodes(channel->getNode0()->getId(), visitedNodes, visitedChannels);
                } else {
                    visitNodes(channel->getNode1()->getId(), visitedNodes, visitedChannels);
                }
            }
        }
    }
}

void Chip::setName(std::string name) {
    this->name = std::move(name);
}

std::string Chip::getName() const {
    return name;
}

int Chip::addChannel(int node0Id, int node1Id, double height, double width, double length, ChannelType type) {
    // create channel
    auto node0 = getOrAddNode(node0Id);
    auto node1 = getOrAddNode(node1Id);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size() + organs.size();
    ;
    auto channel = std::make_unique<Channel>(id, node0, node1, height, width, length, type);

    // add to network as long as channel is still a valid pointer
    network.at(node0->getId()).push_back(channel.get());
    network.at(node1->getId()).push_back(channel.get());

    // add channel
    channels.insert_or_assign(id, std::move(channel));

    return id;
}

int Chip::addChannel(int node0Id, int node1Id, double resistance, ChannelType type) {
    // create channel
    auto node0 = getOrAddNode(node0Id);
    auto node1 = getOrAddNode(node1Id);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size() + organs.size();
    ;
    auto channel = std::make_unique<Channel>(id, node0, node1, resistance, type);

    // add to network as long as channel is still a valid pointer
    network.at(node0->getId()).push_back(channel.get());
    network.at(node1->getId()).push_back(channel.get());

    // add channel
    channels.insert_or_assign(id, std::move(channel));

    return id;
}

int Chip::addMembraneToChannel(int channelId, double height, double width, double poreRadius, double porosity) {
    auto channel = getChannel(channelId);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size() + organs.size();
    auto membrane = std::make_unique<Membrane>(id, channel->getNode0(), channel->getNode1(), height, width, channel->getLength(), poreRadius, porosity);
    membrane->setChannel(channel);

    membranes.insert_or_assign(id, std::move(membrane));

    return id;
}

int Chip::addOrganToMembrane(int membraneId, double height, double width) {
    auto membrane = getMembrane(membraneId);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size() + organs.size();
    auto organ = std::make_unique<Organ>(id, membrane->getNode0(), membrane->getNode1(), height, width, membrane->getLength());
    membrane->setOrgan(organ.get());

    organs.insert_or_assign(id, std::move(organ));

    return id;
}

int Chip::addFlowRatePump(int node0Id, int node1Id, double flowRate) {
    // create pump
    auto node0 = getOrAddNode(node0Id);
    auto node1 = getOrAddNode(node1Id);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size();
    ;
    auto pump = std::make_unique<FlowRatePump>(id, node0, node1, flowRate);

    // add pump
    flowRatePumps.insert_or_assign(id, std::move(pump));

    return id;
}

int Chip::addPressurePump(int node0Id, int node1Id, double pressure) {
    // create pump
    auto node0 = getOrAddNode(node0Id);
    auto node1 = getOrAddNode(node1Id);
    auto id = channels.size() + flowRatePumps.size() + pressurePumps.size() + membranes.size();
    ;
    auto pump = std::make_unique<PressurePump>(id, node0, node1, pressure);

    // add pump
    pressurePumps.insert_or_assign(id, std::move(pump));

    return id;
}

void Chip::addSink(int nodeId) {
    auto sink = getOrAddNode(nodeId);

    //insert sink into sinks (does nothing if sink is already present in sinks)
    sinks.insert(sink);
}

void Chip::addGround(int nodeId) {
    groundNode = getOrAddNode(nodeId);
}

bool Chip::isSink(int nodeId) const {
    return sinks.count(nodes.at(nodeId).get()) == 1;
}

int Chip::getGroundId() const {
    if (groundNode == nullptr) {
        throw std::invalid_argument("Ground node not defined.");
    }
    return groundNode->getId();
}

Node* Chip::getGroundNode() const {
    return groundNode;
}

bool Chip::hasNode(int nodeId) const {
    return nodes.count(nodeId);
}

Node* Chip::getNode(int nodeId) const {
    try {
        return nodes.at(nodeId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Node with ID " + std::to_string(nodeId) + " does not exist.");
    }
}

Channel* Chip::getChannel(int channelId) const {
    try {
        return channels.at(channelId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Channel with ID " + std::to_string(channelId) + " does not exist.");
    }
}

FlowRatePump* Chip::getFlowRatePump(int flowRatePumpId) const {
    try {
        return flowRatePumps.at(flowRatePumpId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Flow rate pump with ID " + std::to_string(flowRatePumpId) + " does not exist.");
    }
}

PressurePump* Chip::getPressurePump(int pressurePumpId) const {
    try {
        return pressurePumps.at(pressurePumpId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Pressure pump with ID " + std::to_string(pressurePumpId) + " does not exist.");
    }
}

Pump* Chip::getPump(int pumpId) const {
    if (pressurePumps.count(pumpId) != 0)
        return getPressurePump(pumpId);
    else if (flowRatePumps.count(pumpId) != 0)
        return getFlowRatePump(pumpId);
    else {
        throw std::invalid_argument("Pump with ID " + std::to_string(pumpId) + " does not exist.");
    }
}

Membrane* Chip::getMembrane(int membraneId) {
    try {
        return membranes.at(membraneId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Membrane with ID " + std::to_string(membraneId) + " does not exist.");
    }
}

Organ* Chip::getOrgan(int organId) {
    try {
        return organs.at(organId).get();
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Organ with ID " + std::to_string(organId) + " does not exist.");
    }
}

const std::unordered_map<int, std::unique_ptr<Channel>>& Chip::getChannels() const {
    return channels;
}

const std::unordered_map<int, std::unique_ptr<Membrane>>& Chip::getMembranes() const {
    return membranes;
}

const std::unordered_map<int, std::unique_ptr<Organ>>& Chip::getOrgans() const {
    return organs;
}

const std::unordered_map<int, std::unique_ptr<Node>>& Chip::getNodes() const {
    return nodes;
}

const std::unordered_map<int, std::unique_ptr<FlowRatePump>>& Chip::getFlowRatePumps() const {
    return flowRatePumps;
}

const std::unordered_map<int, std::unique_ptr<PressurePump>>& Chip::getPressurePumps() const {
    return pressurePumps;
}

const std::vector<Channel*>& Chip::getChannelsAtNode(int nodeId) const {
    try {
        return network.at(nodeId);
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Node with ID " + std::to_string(nodeId) + " does not exist.");
    }
}

Channel* Chip::getChannelBetweenNodes(int nodeId0, int nodeId1) {
    try {
        std::vector<Channel*>& channels = network.at(nodeId0);
        for (auto& channel : channels) {
            if (channel->getNode1()->getId() == nodeId1) {
                return channel;
            }
        }
        throw std::invalid_argument("Channel between ID " + std::to_string(nodeId0) + " and ID " + std::to_string(nodeId1) + " does not exist.");
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("Node with ID " + std::to_string(nodeId0) + " does not exist.");
    }
}

Membrane* Chip::getMembraneBetweenNodes(int nodeId0, int nodeId1) {
    for (auto& [key, membrane] : membranes) {
        if (((membrane->getNode0()->getId() == nodeId0) && (membrane->getNode1()->getId() == nodeId1)) || ((membrane->getNode0()->getId() == nodeId1) && (membrane->getNode1()->getId() == nodeId0))) {
            return membrane.get();
        }
    }
    throw std::invalid_argument("Membrane between ID " + std::to_string(nodeId0) + " and ID " + std::to_string(nodeId1) + " does not exist.");
}

std::vector<Membrane*> Chip::getMembranesAtNode(int nodeId) {
    std::vector<Membrane*> membrane_vector;
    for (auto& [key, membrane] : membranes) {
        if ((membrane->getNode0()->getId() == nodeId) || (membrane->getNode1()->getId() == nodeId)) {
            membrane_vector.push_back(membrane.get());
        }
    }
    return membrane_vector;
}

Organ* Chip::getOrganBetweenNodes(int nodeId0, int nodeId1) {
    for (auto& [key, organ] : organs) {
        if (((organ->getNode0()->getId() == nodeId0) && (organ->getNode1()->getId() == nodeId1)) || ((organ->getNode0()->getId() == nodeId1) && (organ->getNode1()->getId() == nodeId0))) {
            return organ.get();
        }
    }
    throw std::invalid_argument("Organ between ID " + std::to_string(nodeId0) + " and ID " + std::to_string(nodeId1) + " does not exist.");
}

bool Chip::isNetworkValid() {
    // checks if all nodes and channels are connected to ground (if channel network is one graph)
    std::unordered_map<int, bool> visitedNodes;
    std::unordered_map<int, bool> visitedChannels;

    if (nodes.size() == 0) {
        throw std::invalid_argument("No nodes in network.");
    }

    for (auto const& [k, v] : channels) {
        if (v->getLength() <= 0) {
            throw std::invalid_argument("Channel " + std::to_string(k) + ": length is <= 0.");
        }
        if (v->getHeight() <= 0) {
            throw std::invalid_argument("Channel " + std::to_string(k) + ": height is <= 0.");
        }
        if (v->getWidth() <= 0) {
            throw std::invalid_argument("Channel " + std::to_string(k) + ": width is <= 0.");
        }
    }

    for (auto const& [k, v] : nodes) {
        visitedNodes[k] = false;
    }
    for (auto const& [k, v] : channels) {
        visitedChannels[k] = false;
    }

    visitNodes(-1, visitedNodes, visitedChannels);

    std::string errorNodes = "";
    for (auto const& [k, v] : nodes) {
        if (visitedNodes[k] == false) {
            errorNodes.append(" " + std::to_string(k));
        }
    }
    std::string errorChannels = "";
    for (auto const& [k, v] : channels) {
        if (visitedChannels[k] == false) {
            errorChannels.append(" " + std::to_string(k));
        }
    }

    if (errorNodes.length() != 0 || errorChannels.length() != 0) {
        throw std::invalid_argument("Chip is invalid. The following nodes are not connected to ground: " + errorNodes + ". The following channels are not connected to ground: " + errorChannels);
        return false;
    }

    return true;
}

}  // namespace arch