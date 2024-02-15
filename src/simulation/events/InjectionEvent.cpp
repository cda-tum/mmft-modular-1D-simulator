#include "InjectionEvent.h"

#include "DropletBoundary.h"

namespace sim {

InjectionEvent::InjectionEvent(double time, Injection& injection) : Event(time, 1), injection(injection) {}

void InjectionEvent::performEvent() {
    // injection position of the droplet (center of the droplet)
    auto channelPosition = injection.getInjectionPosition();
    auto channel = channelPosition.getChannel();
    auto droplet = injection.getDroplet();

    // for the injection the two boundaries (basically a head and a tail) of the droplet must lie inside the channel (the volume of the droplet must be small enough)
    // this is already checked at the creation of the injection and, thus, we don't need this check here

    // the droplet length is a relative value between 0 and 1
    double dropletLength = droplet->getVolume() / channel->getVolume();

    // compute position of the two boundaries
    double position0 = (channelPosition.getPosition() - dropletLength / 2);  // is always the position which lies closer to 0
    double position1 = (channelPosition.getPosition() + dropletLength / 2);  // is always the position which lies closer to 1

    // create corresponding boundaries
    // since boundary0 always lies closer to 0, the volume of this boundary points to node1
    // since boundary1 always lies closer to 1, the volume of this boundary points to node0
    droplet->addBoundary(channel, position0, false, BoundaryState::NORMAL);
    droplet->addBoundary(channel, position1, true, BoundaryState::NORMAL);

    // set droplet state
    droplet->setDropletState(DropletState::NETWORK);
}

}  // namespace sim