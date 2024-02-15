#include "ResistanceModels.h"

#include <math.h>

namespace sim {

// ### ResistanceModel0 ###
ResistanceModel0::ResistanceModel0(double continuousPhaseViscosity) : continuousPhaseViscosity(continuousPhaseViscosity) {}

double ResistanceModel0::getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const {
    double a = computeFactorA(channel->getWidth(), channel->getHeight());

    return channel->getLength() * a * continuousPhaseViscosity / (channel->getWidth() * pow(channel->getHeight(), 3));
}

double ResistanceModel0::computeFactorA(double width, double height) const {
    return 12 / (1 - 192 * height * tanh(M_PI * width / (2 * height)) / (pow(M_PI, 5) * width));
}

double ResistanceModel0::getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const {
    double a = computeFactorA(channel->getWidth(), channel->getHeight());
    double dropletLength = volumeInsideChannel / (channel->getWidth() * channel->getHeight());

    return 3 * dropletLength * a * continuousPhaseViscosity / (channel->getWidth() * pow(channel->getHeight(), 3));
}

// ### ResistanceModel1 ###
ResistanceModel1::ResistanceModel1() {}

double ResistanceModel1::getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const {
    return channel->getVolume();
}

double ResistanceModel1::getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const {
    return 3 * channel->getWidth() * channel->getHeight() * volumeInsideChannel;
}

// ### ResistanceModel2 ###
ResistanceModel2::ResistanceModel2(double continuousPhaseViscosity) : continuousPhaseViscosity(continuousPhaseViscosity) {}

double ResistanceModel2::getChannelResistance(arch::Channel const* const channel, double avgViscosityInChannel) const {
    double a = computeFactorA(channel->getWidth(), channel->getHeight());

    return channel->getLength() * a * avgViscosityInChannel / (channel->getWidth() * pow(channel->getHeight(), 3));
}

double ResistanceModel2::computeFactorA(double width, double height) const {
    return 12 / (1 - 192 * height * tanh(M_PI * width / (2 * height)) / (pow(M_PI, 5) * width));
}

double ResistanceModel2::getDropletResistance(arch::Channel const* const channel, sim::Droplet* droplet, double volumeInsideChannel, double avgViscosityInChannel) const {
    double a = computeFactorA(channel->getWidth(), channel->getHeight());
    double dropletLength = volumeInsideChannel / (channel->getWidth() * channel->getHeight());

    return 3 * dropletLength * a * avgViscosityInChannel / (channel->getWidth() * pow(channel->getHeight(), 3));
}

}  // namespace sim