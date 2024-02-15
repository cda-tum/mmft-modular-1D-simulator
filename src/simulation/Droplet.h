/**
 * @file Droplet.h
 */

#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Channel.h"
#include "DropletBoundary.h"
#include "Fluid.h"
#include "IResistanceModel.h"

namespace sim {

/**
 * @brief Enum to specify in which state the droplet currently is in.
 */
enum class DropletState {
    INJECTION,  ///< Droplet planned to be injected but currently not yet in the network.
    NETWORK,    ///< Droplet currently flows through the network.
    TRAPPED,    ///< Droplet is trapped in the network.
    SINK        ///< Droplet has left the network (is in the sink).
};

/**
 * @brief Class to specify a droplet.
 */
class Droplet {
  private:
    int const id;                                              ///< Unique identifier of the droplet.
    std::string name = "";                                     ///< Name of the droplet.
    double volume;                                             ///< Volume of the droplet in m^3.
    Fluid* fluid;                                              ///< Pointer to fluid of which the droplet consists of.
    std::vector<Droplet*> mergedDroplets;                      ///< List of previous droplets, if this droplet got merged.
    DropletState dropletState = DropletState::INJECTION;       ///< Current state of the droplet
    std::vector<std::unique_ptr<DropletBoundary>> boundaries;  ///< Boundaries of the droplet
    std::vector<arch::Channel*> channels;                      ///< Contains the channels, that are completely occupied by the droplet (can happen in short channels or with large droplets).

  public:
    /**
     * @brief Specify a droplet.
     * @param[in] id Unique identifier of the droplet.
     * @param[in] volume Volume of the droplet in m^3.
     * @param[in] fluid Pointer to fluid the droplet consists of.
     */
    Droplet(int id, double volume, Fluid* fluid);

    /**
     * @brief Change volume of droplet.
     * @param[in] volume New volume of the droplet.
     */
    void setVolume(double volume);

    /**
     * @brief Set name of droplet.
     * @param[in] name New name of droplet.
     */
    void setName(std::string name);

    /**
     * @brief Change state in which the droplet currently is in.
     * @param[in] dropletState The new state in which the droplet is in.
     */
    void setDropletState(DropletState dropletState);

    /**
     * @brief Returns unique identifier of the droplet.
     * @return Unique identifier of the droplet.
     */
    int getId() const;

    /**
     * @brief Retrieve the name of the droplet.
     * @return The name of the droplet.
     */
    std::string getName() const;

    /**
     * @brief Retrieve the volume of droplet.
     * @return The volume of the droplet in m^3.
     */
    double getVolume() const;

    /**
     * @brief Retrieve the state in which the droplet is currently in.
     * @return The droplet state the droplet has at the moment.
     */
    DropletState getDropletState() const;

    /**
     * @brief Retrieve the fluid of the droplet.
     * @return The fluid the droplet consists of in m^3.
     */
    const Fluid* getFluid() const;

    /**
     * @brief Add the resistance the droplet causes to the channels it currently occupies.
     * @param[in] model The resistance model on which basis the resistance caused by the droplet is calculated.
     */
    void addResistances(const IResistanceModel& model);

    /**
     * @brief Get the Boundaries object
     * @return all boundaries
     */
    const std::vector<std::unique_ptr<sim::DropletBoundary>>& getBoundaries() const;

    /**
     * @brief Get all fully occupied channels
     * @return all fully occupied channels
     */
    std::vector<arch::Channel*>& getFullyOccupiedChannels();

    /**
     * @brief If the droplet currently is at a bifurcation
     * @return true
     * @return false
     */
    bool isAtBifurcation();

    /**
     * @brief If the droplet currently is inside a single channel
     * @return true
     * @return false
     */
    bool isInsideSingleChannel();

    /**
     * @brief Add a boundary to the boundary list of the droplet
     * @param channel Channel
     * @param position Position within channel
     * @param volumeTowardsNode0 Direction in which the droplet lies within the channel (in regards to node0)
     * @param state State the boundary is in
     */
    void addBoundary(arch::Channel* channel, double position, bool volumeTowardsNode0, BoundaryState state);

    /**
     * @brief Add fully occupied channel to the fully occupied channel list.
     * @param channel New fully occupied channel.
     */
    void addFullyOccupiedChannel(arch::Channel* channel);

    /**
     * @brief Remove boundary from the boundary list.
     * @param boundaryReference Reference to the boundary that should be removed.
     */
    void removeBoundary(DropletBoundary& boundaryReference);

    /**
     * @brief Remove fully occupied channel from the fully occupied channel list.
     * @param channelId Id of the channel that should be removed.
     */
    void removeFullyOccupiedChannel(int channelId);

    /**
     * @brief Get a list of boundaries that are "connected" to the corresponding node
     * @param nodeId Id of the node
     * @param doNotConsider boundary that should not be included in the list (helpful, during request for droplet merging)
     * @return List of connected boundaries
     */
    std::vector<DropletBoundary*> getConnectedBoundaries(int nodeId, DropletBoundary* doNotConsider = nullptr);

    /**
     * @brief Get a list of fully occupied channels that are "connected" to the corresponding node 
     * @param nodeId Id of the node
     * @return List of connected fully occupied channels
     */
    std::vector<arch::Channel*> getConnectedFullyOccupiedChannels(int nodeId);

    /**
     * @brief Update the flow-rates of the droplet boundaries according to the flowRates inside the channels
     * @param chip Chip
     * @param slipFactor indicates the ratio between the velocity of the droplet boundary and the continuous phase (droplet is faster than the continuous phase)
     */
    void updateBoundaries(const arch::Chip& chip, double slipFactor);

    /**
     * @brief Adds a droplet from which the current droplet was created by merging.
     * @param droplet Pointer to droplet.
     */
    void addMergedDroplet(Droplet* droplet);

    /**
     * @brief Gets the list of droplets, that created the actual droplet due to merging.
     * @return List of merged droplets.
     */
    const std::vector<Droplet*>& getMergedDroplets() const;
};

}  // namespace sim