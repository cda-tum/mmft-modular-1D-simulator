/**
 * @file Chip.h
 */
#pragma once

#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "Channel.h"
#include "FlowRatePump.h"
#include "Membrane.h"
#include "Node.h"
#include "Organ.h"
#include "PressurePump.h"

namespace arch {

/**
 * @brief Class to specify a microfluidic chip with all of its components.
 */
class Chip {
  private:
    std::string name;                                                      ///< Name of the chip.
    std::unordered_map<int, std::unique_ptr<Node>> nodes;                  ///< Nodes the network of the chip consists of.
    std::set<Node*> sinks;                                                 ///< Ids of nodes that are sinks.
    Node* groundNode = nullptr;                                            ///< Pointer to the ground node.
    std::unordered_map<int, std::unique_ptr<Channel>> channels;            ///< Map of ids and channel pointer of all channels the chip consists of.
    std::unordered_map<int, std::unique_ptr<Membrane>> membranes;          ///< Map of ids and membrane pointer of all membranes the chip consists of.
    std::unordered_map<int, std::unique_ptr<Organ>> organs;                ///< Map of ids and organ pointer of all organs the chip consists of.
    std::unordered_map<int, std::unique_ptr<PressurePump>> pressurePumps;  ///< Map of ids and pressure pump pointers of all pressure pumps the chip consists of.
    std::unordered_map<int, std::unique_ptr<FlowRatePump>> flowRatePumps;  ///< Map of ids of flow rate pump pointers of all flow rate pumps the chip consists of.
    std::unordered_map<int, std::vector<Channel*>> network;                ///< Network of nodes and corresponding channels at these nodes on the chip.

    /**
     * @brief To add components (vertices) to the chip, the nodes of these vertices need to be specified. This function tries to find the specified node or creates a new one with this id if it does not exist yet and returns it.
     * @param[in] nodeId Id of the node to be returned.
     * @return Node with the given id.
     */
    Node* getOrAddNode(int nodeId);

    /**
     * @brief Goes through network and sets all nodes and channels that are visited to true.
     * @param[in] id Id of the node that is visited.
     * @param[in, out] visitedNodes Reference to a map that stores which nodes have already been visited.
     * @param[in, out] visitedChannels Reference to a map that stores which channels have already been visited.
     */
    void visitNodes(int id, std::unordered_map<int, bool>& visitedNodes, std::unordered_map<int, bool>& visitedChannels);

  public:
    /**
     * @brief Constructor of the chip
     */
    Chip();

    /**
     * @brief Sets the name of the chip.
     * @param[in] name Name of the chip.
     */
    void setName(std::string name);

    /**
     * @brief Returns the name of the chip.
     * @return Name of the chip.
     */
    std::string getName() const;

    /**
     * @brief Adds a new channel to the chip.
     * @param[in] node0Id Id of the node at one end of the channel.
     * @param[in] node1Id Id of the node at the other end of the channel.
     * @param[in] height Height of the channel in m.
     * @param[in] width Width of the channel in m.
     * @param[in] length Length of the channel in m.
     * @param[in] type What kind of channel it is.
     * @return Id of the newly created channel.
     */
    int addChannel(int node0Id, int node1Id, double height, double width, double length, ChannelType type);

    /**
     * @brief Adds a new channel to the chip.
     * @param[in] node0Id Id of the node at one end of the channel.
     * @param[in] node1Id Id of the node at the other end of the channel.
     * @param[in] resistance Resistance of the channel in Pas/m^3.
     * @param[in] type What kind of channel it is.
     * @return Id of the newly created channel.
     */
    int addChannel(int node0Id, int node1Id, double resistance, ChannelType type);

    /**
     * @brief Creates and adds a membrane to a channel in the simulator.
     * @param[in] channelId Id of the channel. Channel defines nodes, length and width.
     * @param[in] height Height of the channel in m.
     * @param[in] poreSize Size of the pores in m.
     * @param[in] porosity Porosity of the membrane in % (between 0 and 1).
     * @return Id of the membrane.
     */
    int addMembraneToChannel(int channelId, double height, double width, double poreRadius, double porosity);

    /**
     * @brief Creates and adds a organ to a membrane in the simulator.
     * @param[in] membraneId Id of the membrane. Membrane defines nodes, length and width.
     * @param[in] height Height of the organ in m.
     */
    int addOrganToMembrane(int membraneId, double height, double width);

    /**
     * @brief Adds a new flow rate pump to the chip.
     * @param[in] node0Id Id of the node at one end of the flow rate pump.
     * @param[in] node1Id Id of the node at the other end of the flow rate pump.
     * @param[in] flowRate Volumetric flow rate of the pump in m^3/s.
     * @return Id of the newly created flow rate pump.
     */
    int addFlowRatePump(int node0Id, int node1Id, double flowRate);

    /**
     * @brief Adds a new pressure pump to the chip.
     * @param[in] node0Id Id of the node at one end of the pressure pump.
     * @param[in] node1Id Id of the node at the other end of the pressure pump.
     * @param[in] pressure Pressure of the pump in Pas/m^3.
     * @return Id of the newly created pressure pump.
     */
    int addPressurePump(int node0Id, int node1Id, double pressure);

    /**
     * @brief Specifies a node as sink.
     * @param[in] nodeId Id of the node that is a sink.
     */
    void addSink(int nodeId);

    /**
     * @brief Adds or sets a node as the ground node, i.e., this node has a pressure value of 0 and acts as a reference node for all other nodes.
     * @param[in] nodeId Id of the node that should be the ground node of the network.
     */
    void addGround(int nodeId);

    /**
     * @brief Checks and returns if a node is a sink.
     * @param[in] nodeId Id of the node that should be checked.
     * @return If the node with the specified id is a sink.
     */
    bool isSink(int nodeId) const;

    /**
     * @brief Returns the id of the ground node.
     * @return Id of the ground node.
     */
    int getGroundId() const;

    /**
     * @brief Returns a pointer to the ground node.
     * @return Pointer to the ground node.
     */
    Node* getGroundNode() const;

    /**
     * @brief Checks if a node with the specified id exists in the network.
     * @param[in] nodeId Id of the node to check.
     * @return If such a node exists.
     */
    bool hasNode(int nodeId) const;

    /**
     * @brief Get pointer to node with the specified id.
     * @param[in] nodeId Id of the node to get.
     * @returns Pointer to the node with this id.
     */
    Node* getNode(int nodeId) const;

    /**
     * @brief Get pointer to channel with the specified id.
     * @param[in] channelId Id of the channel to get.
     * @return Pointer to the channel with this id.
     */
    Channel* getChannel(int channelId) const;

    /**
     * @brief Get pointer to flow rate pump with the specified id.
     * @param[in] flowRatePumpId Id of the flow rate pump to get.
     * @return Pointer to the flow rate pump with this id.
     */
    FlowRatePump* getFlowRatePump(int flowRatePumpId) const;

    /**
     * @brief Get pointer to pressure pump with the specified id.
     * @param[in] pressurePumpId Id of the pressure pump to get.
     * @return Pointer to the pressure pump with this id.
     */
    PressurePump* getPressurePump(int pressurePumpId) const;

    /**
     * @brief Get pointer to a pump with the specified id.
     * 
     * @param pumpId Id of the pump.
     * @return Pointer to the pump with this id.
     */
    Pump* getPump(int pumpId) const;

    /**
     * @brief Get pointer to a membrane with the specified id.
     * 
     * @param membraneId Id of the membrane.
     * @return Pointer to the membrane with this id.
     */
    Membrane* getMembrane(int membraneId);

    /**
     * @brief Get pointer to an organ with the specified id.
     * 
     * @param organId Id of the organ.
     * @return Pointer to the organ with this id.
     */
    Organ* getOrgan(int organId);

    /**
     * @brief Get a map of all channels of the chip.
     @return Map that consists of the channel ids and pointers to the corresponding channels.
     */
    const std::unordered_map<int, std::unique_ptr<Channel>>& getChannels() const;

    /**
     * @brief Get a map of all membranes of the chip.
     @return Map that consists of the membrane ids and pointers to the corresponding membranes.
     */
    const std::unordered_map<int, std::unique_ptr<Membrane>>& getMembranes() const;

    /**
     * @brief Get a map of all organs of the chip.
     @return Map that consists of the organ ids and pointers to the corresponding organs.
     */
    const std::unordered_map<int, std::unique_ptr<Organ>>& getOrgans() const;

    /**
     * @brief Get a map of all nodes of the chip.
     * @return Map that consists of the node ids and pointers to the corresponding nodes.
     */
    const std::unordered_map<int, std::unique_ptr<Node>>& getNodes() const;

    /**
     * @brief Get a map of all flow rate pumps of the chip.
     * @return Map that consists of the flow rate pump ids and pointers to the corresponding flow rate pumps.
     */
    const std::unordered_map<int, std::unique_ptr<FlowRatePump>>& getFlowRatePumps() const;

    /**
     * @brief Get a map of all pressure pumps of the chip.
     * @return Map that consists of the pressure pump ids and pointers to the corresponding pressure pumps.
     */
    const std::unordered_map<int, std::unique_ptr<PressurePump>>& getPressurePumps() const;

    /**
     * @brief Get a map of all channels at a specific node.
     * @param[in] nodeId Id of the node at which the adherent channels should be returned.
     * @return Vector of pointers to channels adherent to this node.
     */
    const std::vector<Channel*>& getChannelsAtNode(int nodeId) const;

    /**
     * @brief Get all edges at a specific node.
     * 
     * @param nodeId Id of the node.
     * @return Vector of pointers to all edges connected to this node.
     */
    const std::vector<Edge*>& getEdgesAtNode(int nodeId);

    /**
     * @brief Get the channel that are connected to both specified nodes.
     * 
     * @param nodeId0 Id of node 0.
     * @param nodeId1 Id of node 1.
     * @return Pointer to the channel lies between these nodes.
     */
    Channel* getChannelBetweenNodes(int nodeId0, int nodeId1);

    /**
     * @brief Get the membrane that is connected to both specified nodes.
     * 
     * @param nodeId0 Id of node 0.
     * @param nodeId1 Id of node 1.
     * @return Pointer to the membrane that lies between these nodes.
     */
    Membrane* getMembraneBetweenNodes(int nodeId0, int nodeId1);

    /**
     * @brief Get vector of all membranes that are connected to the specified node.
     * 
     * @param nodeId Id of the node.
     * @return Vector containing pointers to all membranes that are connected to this node.
     */
    std::vector<Membrane*> getMembranesAtNode(int nodeId);

    /**
     * @brief Get the organ that lies between two nodes.
     * 
     * @param nodeId0 Id of node0.
     * @param nodeId1 Id of node1.
     * @return Pointer to the organ that lies between the two nodes.
     */
    Organ* getOrganBetweenNodes(int nodeId0, int nodeId1);

    /**
     * @brief Checks if chip network is valid in the sense that all nodes and channels need to be connected to ground (and channel network must be one graph).
     * @return If the network is valid.
     */
    bool isNetworkValid();
};

}  // namespace arch