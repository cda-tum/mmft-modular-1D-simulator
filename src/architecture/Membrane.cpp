#include "Membrane.h"

#include "Edge.h"
#include "Node.h"
#include "ResistanceModels.h"

#include <cmath>
#include <stdexcept>

namespace arch {

Membrane::Membrane(int id, Node* node0, Node* node1, double height, double width, double length, double poreRadius, double porosity) : Edge(id, node0, node1), height(height), width(width), length(length), poreRadius(poreRadius), porosity(porosity) {}

Membrane::Membrane(int id, Node* node0, Node* node1, double resistance) : Edge(id, node0, node1), membraneResistance(resistance) {}

void Membrane::setDimensions(double width, double height, double length) {
    this->height = height;
    this->width = width;
    this->length = length;
}

void Membrane::setHeight(double height) {
    this->height = height;
}

void Membrane::setWidth(double width) {
    this->width = width;
}

void Membrane::setLength(double length) {
    this->length = length;
}

void Membrane::setPoreRadius(double poreRadius) {
    this->poreRadius = poreRadius;
}

void Membrane::setPorosity(double porosity) {
    this->porosity = porosity;
}

void Membrane::setChannel(Channel* channel) {
    this->channel = channel;
}

void Membrane::setOrgan(Organ* organ) {
    this->organ = organ;
}

double Membrane::getHeight() const {
    return height;
}

double Membrane::getWidth() const {
    return width;
}

double Membrane::getLength() const {
    return length;
}

Channel* Membrane::getChannel() const {
    return channel;
}

Organ* Membrane::getOrgan() const {
    return organ;
}

double Membrane::getPoreRadius() const {
    return poreRadius;
}

double Membrane::getPoreDiameter() const {
    return poreRadius * 2;
}

double Membrane::getPorosity() const {
    return porosity;
}

double Membrane::getNumberOfPores(double area) const {
    return (porosity * area) / (M_PI * pow(poreRadius, 2));
}

double Membrane::getPoreDensity() const {
    return (porosity * M_PI * pow(poreRadius, 2));
}

double Membrane::getArea() const {
    return width * length;
}

double Membrane::getVolume() const {
    return width * height * length;
}

double Membrane::getPressure() const {
    return node0->getPressure() - node1->getPressure();
}

double Membrane::getFlowRate() const {
    return getPressure() / getResistance();
}

double Membrane::getResistance() const {
    return membraneResistance;
}

double Membrane::getConcentrationChange(double resistance, double timeStep, double concentrationDifference, double currTime) const {
#define f(time, concentration, permeability) (permeability * concentration)
    auto* membraneChannel = this->getChannel();
    double permeability = 1 / resistance;
    double k1 = timeStep * (f(currTime, concentrationDifference, permeability));
    double k2 = timeStep * (f((currTime + timeStep / 2), (concentrationDifference + k1 / 2), permeability));
    double k3 = timeStep * (f((currTime + timeStep / 2), (concentrationDifference + k2 / 2), permeability));
    double k4 = timeStep * (f((currTime + timeStep), (concentrationDifference + k3), permeability));
    double concentrationChange = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return concentrationChange;
}

}  // namespace arch