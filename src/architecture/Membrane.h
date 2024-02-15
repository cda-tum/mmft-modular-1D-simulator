/**
 * @file Membrane.h
 */

#pragma once

#include "Channel.h"
#include "Edge.h"
#include "IResistance.h"
#include "Node.h"
#include "Organ.h"

namespace arch {

/**
 * @brief Class to specify a membrane, which is a component of a chip in which droplet can flow.
 */
class Membrane : public virtual Edge, public virtual nodal::IResistance {
  private:
    double height = 0;              ///< Height of a membrane in m.
    double width = 0;               ///< Width of a membrane in m.
    double length = 0;              ///< Length of a membrane in m.
    double poreRadius = 0;          ///< Radius of the pores.
    double porosity = 0;            ///< Porosity of the membrane in % (between 0.0 and 1.0)
    double numberOfPores = 0;       ///< Numbers of pores of the membrane.
    double membraneResistance = 0;  ///< Resistance of the membrane in Pas/m^3.
    Channel* channel;               ///< Membrane on which the barrier is attached (length must be equal).
    Organ* organ;                   ///< Pointer to the organ this membrane is attached to (length must be equal).

  public:
    /**
     * @brief Constructor of a membrane.
     * @param[in] id Id of the membrane.
     * @param[in] node0 Node at one end of the membrane.
     * @param[in] node1 Node at other end of the membrane.
     * @param[in] height Height of the membrane in m.
     * @param[in] width Width of the membrane in m.
     * @param[in] length Length of the membrane in m.
     * @param[in] poreRadius Radius of the pores in the membrane in m.
     * @param[in] porosity Porosity of the membrane.
     * @param[in] type Type of the membrane.
     */
    Membrane(int id, Node* node0, Node* node1, double height, double width, double length, double poreRadius, double porosity);

    /**
     * @brief Constructor of a membrane.
     * @param[in] id Id of the membrane.
     * @param[in] node0 Node at one end of the membrane.
     * @param[in] node1 Node at other end of the membrane.
     * @param[in] resistance Resistance of the membrane in Pas/m^3.
     * @param[in] type Type of the membrane.
     */
    Membrane(int id, Node* node0, Node* node1, double resistance);

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
     * @brief Set radius of the pores of the membrane.
     * 
     * @param poreRadius Radius of the pores in m.
     */
    void setPoreRadius(double poreRadius);

    /**
     * @brief Set porosity of the membrane.
     * 
     * @param porosity Porosity in % (between 0.0 and 1.0)
     */
    void setPorosity(double porosity);

    /**
     * @brief Set channel the membrane is connected to.
     * 
     * @param channel Pointer to the channel the membrane is attached to.
     */
    void setChannel(Channel* channel);

    /**
     * @brief Set organ the membrane is connected to.
     * 
     * @param organ Pointer to the organ the membrane is attached to.
     */
    void setOrgan(Organ* organ);

    /**
     * @brief Get the concentration change caused by the membrane for a specific species in a mixture given the concentration difference between the channel and the organ tank to which the membrane is connected to.
     * 
     * @param resistance Resistance of the membrane for the mixture area for a specific species.
     * @param timeStep Time step for which the simulation is forwarded.
     * @param concentrationDifference Concentration difference between the channel and the organ tank.
     * @param currTime Current simulation time.
     * @return Absolute concentration change in mol.
     */
    double getConcentrationChange(double resistance, double timeStep, double concentrationDifference, double currTime) const;

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
     * @brief Get pointer to the channel the membrane is connected to.
     * 
     * @return Pointer to the channel the membrane is connected to.
     */
    Channel* getChannel() const;

    /**
    * @brief Get pointer to the organ the membrane is connected to.
    * 
    * @return Pointer to the organ the membrane is connected to.
    */
    Organ* getOrgan() const;

    /**
     * @brief Get radius of the pores of the membrane.
     * 
     * @return Radius of a pore of the membrane in m. 
     */
    double getPoreRadius() const;

    /**
     * @brief Get diameter of the pores of the membrane.
     * 
     * @return Diameter of a pore of the membrane in m.
     */
    double getPoreDiameter() const;

    /**
     * @brief Get the porosity of the membrane.
     * 
     * @return Porosity of the membrane in % (between 0.0 and 1.0).
     */
    double getPorosity() const;

    /**
     * @brief Get the number of pores of the membrane.
     * @param area Area for which the number of pores should be calculated. This is used to get the number of pores for a specific area covered by one mixture.
     * @return The number of pores of the membrane.
     */
    double getNumberOfPores(double area) const;

    /**
     * @brief Get the density of the pores on the membrane.
     * 
     * @return Density of the pores.
     */
    double getPoreDensity() const;

    /**
     * @brief Returns area of a membrane.
     * @returns Area in m^2.
     */
    double getArea() const;

    /**
     * @brief Calculates and returns volume of the membrane.
     * @returns Volume of a channel in m^3.
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