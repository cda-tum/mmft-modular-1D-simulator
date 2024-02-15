#include "Organ.h"

#include "Edge.h"
#include "Node.h"

#include <cmath>
#include <stdexcept>

namespace arch {

Organ::Organ(int id, Node* node0, Node* node1, double height, double width, double length) : Edge(id, node0, node1), height(height), width(width), length(length) {}

Organ::Organ(int id, Node* node0, Node* node1, double resistance) : Edge(id, node0, node1), organResistance(resistance) {}

void Organ::setDimensions(double width, double height, double length) {
    this->height = height;
    this->width = width;
    this->length = length;
}

void Organ::setHeight(double height) {
    this->height = height;
}

void Organ::setWidth(double width) {
    this->width = width;
}

void Organ::setLength(double length) {
    this->length = length;
}

double Organ::getHeight() const {
    return height;
}

double Organ::getWidth() const {
    return width;
}

double Organ::getLength() const {
    return length;
}

double Organ::getArea() const {
    return width * length;
}

double Organ::getVolume() const {
    return width * height * length;
}

double Organ::getPressure() const {
    return node0->getPressure() - node1->getPressure();
}

double Organ::getFlowRate() const {
    // there is no flow in an organ, thus the organ does not have a flow-rate
    // to ensure time-accurate concentration changes
    // mixtures move through the organ based on the flow-rate of the connected channel
    return 0;
}

double Organ::getResistance() const {
    return organResistance;
}

}  // namespace arch