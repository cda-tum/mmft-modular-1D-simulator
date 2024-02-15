
#include "MergeChannelEvent.h"

#include "Droplet.h"
#include "DropletBoundary.h"
#include "Simulation.h"

namespace sim {

MergeChannelEvent::MergeChannelEvent(double time, Droplet& droplet0, DropletBoundary& boundary0, Droplet& droplet1, DropletBoundary& boundary1, Simulation& simulation) : Event(time, 0), droplet0(droplet0), boundary0(boundary0), droplet1(droplet1), boundary1(boundary1), simulation(simulation) {}

void MergeChannelEvent::performEvent() {
    auto newDroplet = simulation.mergeDroplets(droplet0.getId(), droplet1.getId());

    // add boundaries from droplet0
    for (auto& boundary : droplet0.getBoundaries()) {
        // do not add boundary0 to new droplet
        if (boundary.get() == &boundary0) {
            continue;
        }

        // add boundary
        newDroplet->addBoundary(boundary->getChannelPosition().getChannel(), boundary->getChannelPosition().getPosition(), boundary->isVolumeTowardsNode0(), boundary->getState());
    }

    // add fully occupied channels from droplet0
    for (auto& channel : droplet0.getFullyOccupiedChannels()) {
        newDroplet->addFullyOccupiedChannel(channel);
    }

    // add boundaries from droplet1
    for (auto& boundary : droplet1.getBoundaries()) {
        // do not add boundary1 to new droplet
        if (boundary.get() == &boundary1) {
            continue;
        }

        // add boundary
        newDroplet->addBoundary(boundary->getChannelPosition().getChannel(), boundary->getChannelPosition().getPosition(), boundary->isVolumeTowardsNode0(), boundary->getState());
    }

    // add fully occupied channels from droplet1
    for (auto& channel : droplet1.getFullyOccupiedChannels()) {
        newDroplet->addFullyOccupiedChannel(channel);
    }

    // check if droplet0 and droplet1 are not inside a single channel, because then a fully occupied channel has to be added
    if (!droplet0.isInsideSingleChannel() && !droplet1.isInsideSingleChannel()) {
        // boundary0 and boundary1 must have the same channel
        newDroplet->addFullyOccupiedChannel(boundary0.getChannelPosition().getChannel());
    }

    // add/remove droplets from network
    newDroplet->setDropletState(DropletState::NETWORK);
    droplet0.setDropletState(DropletState::SINK);
    droplet1.setDropletState(DropletState::SINK);
}

}  // namespace sim