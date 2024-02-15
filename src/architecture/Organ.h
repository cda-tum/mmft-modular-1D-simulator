/**
 * @file Organ.h
 */

#pragma once

#include "Channel.h"
#include "Edge.h"
#include "IResistance.h"
#include "Node.h"

namespace arch {

/**
 * @brief Class to specify a membrane, which is a component of a chip in which droplet can flow.
 */
class Organ : public virtual Edge, public virtual nodal::IResistance {
  private:
    double height = 0;           ///< Height of a membrane in m.
    double width = 0;            ///< Width of a membrane in m.
    double length = 0;           ///< Length of a membrane in m.
    double organResistance = 0;  ///< Resistance of the membrane in Pas/m^3.

  public:
    /**
     * @brief Constructor of a membrane.
     * @param[in] id Id of the membrane.
     * @param[in] node0 Node at one end of the membrane.
     * @param[in] node1 Node at other end of the membrane.
     * @param[in] height Height of the membrane in m.
     * @param[in] width Width of the membrane in m.
     * @param[in] length Length of the membrane in m.
     */
    Organ(int id, Node* node0, Node* node1, double height, double width, double length);

    /**
     * @brief Constructor of a membrane.
     * @param[in] id Id of the membrane.
     * @param[in] node0 Node at one end of the membrane.
     * @param[in] node1 Node at other end of the membrane.
     * @param[in] resistance Resistance of the membrane in Pas/m^3.
     * @param[in] type Type of the membrane.
     */
    Organ(int id, Node* node0, Node* node1, double resistance);

    /**
     * @brief Set dimensions of a membrane.
     * @param[in] width Width of a membrane in m.
     * @param[in] height Height of a membrane in m.
     * @param[in] length Length of a membrane in m.
     */
    void setDimensions(double width, double height, double length);

    /**
     * @brief Set height of a membrane.
     * @param height New height of membrane in m.
     */
    void setHeight(double height);

    /**
     * @brief Set width of membrane.
     * @param[in] width New width of a membrane in m.
     */
    void setWidth(double width);

    /**
     * @brief Set length of membrane.
     * @param[in] length New length of this membrane in m.
     */
    void setLength(double length);

    /**
     * @brief Returns the height of this membrane.
     * @returns Height of channel in m.
     */
    double getHeight() const;

    /**
     * @brief Returns the width of this membrane.
     * @returns Width of channel in m.
     */
    double getWidth() const;

    /**
     * @brief Returns the length of this membrane.
     * @returns Length of channel in m.
     */
    double getLength() const;

    /**
     * @brief Returns area of a membrane.
     * @returns Area in m^2.
     */
    double getArea() const;

    /**
     * @brief Calculates and returns volume of the membrane.
     * @returns Volumne of a channel in m^3.
     */
    double getVolume() const;

    /**
     * @brief Calculates and returns pressure difference over a channel.
     * @returns Pressure in Pa.
     */
    double getPressure() const override;

    /**
     * @brief Calculate flow rate within the channel.
     * @returns Flow rate in m^3/s.
     */
    double getFlowRate() const override;

    /**
     * @brief Calculate and returns overall resistance caused by the channel itself and the droplets within the channel.
     * @returns Overall resistance in Pas/m^3.
     */
    double getResistance() const override;
};

}  // namespace arch