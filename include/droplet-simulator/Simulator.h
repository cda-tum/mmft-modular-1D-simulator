/**
 * @file Simulator.h
 */

#pragma once

#include <memory>

#include "Results.h"

namespace droplet {

/**
 * @brief Implementation of the simulator.
 */
class SimulatorImpl;

/**
 * @brief Public interface of the simulator.
 */
class Simulator {
  private:
    std::unique_ptr<SimulatorImpl> impl;

  public:
    Simulator();
    ~Simulator();

    /**
     * @brief Creates and adds a normal channel to the simulator.
     * @param[in] node0Id Id of the node at one end of the channel.
     * @param[in] node1Id If of the node at other end of the channel
     * @param[in] height Height of the channel in m.
     * @param[in] width Width of the channel in m.
     * @param[in] length Length of the channel in m.
     * @return Id of the channel.
     */
    int addChannel(int node0Id, int node1Id, double height, double width, double length);

    /**
     * @brief Creates and adds a membrane to a channel in the simulator.
     * @param[in] channelId Id of the channel. Channel defines nodes and length.
     * @param[in] height Height of the membrane in m.
     * @param[in] width Width of the membrane in m.
     * @param[in] poreRadius Radius of the pores in m.
     * @param[in] porosity Porosity of the membrane in % (between 0 and 1).
     * @return Id of the membrane.
     */
    int addMembraneToChannel(int channelId, double height, double width, double poreRadius, double porosity);

    /**
     * @brief Creates and adds an organ to a membrane in the simulator.
     * 
     * @param membraneId Id of the membrane the organ should be added to.
     * @param height Height of the organ in m.
     * @param width Width of the organ in m.
     * @return Id of the organ.
     */
    int addOrganToMembrane(int membraneId, double height, double width);

    /**
     * @brief Creates and adds bypass channel to the simulator.
     * @param[in] node0Id Id of the node at one end of the channel.
     * @param[in] node1Id Id of the node at the other end of the channel.
     * @param[in] height Height of the channel in m.
     * @param[in] width Width of the channel in m.
     * @param[in] length Length of the channel in m.
     * @return Id of the channel.
     */
    int addBypassChannel(int node0Id, int node1Id, double height, double width, double length);

    /**
     * @brief Creates and adds flow rate pump to the simulator.
     * @param[in] node0Id Id of the node at one end of the flow rate pump.
     * @param[in] node1Id Id of the node at the other end of the flow rate pump.
     * @param[in] flowRate Flow rate of the pump in m^3/s.
     * @return Id of the flow rate pump.
     */
    int addFlowRatePump(int node0Id, int node1Id, double flowRate);

    /**
     * @brief Creates and adds pressure pump to the simulator.
     * @param[in] node0Id Id of the node at one end of the pressure rate pump.
     * @param[in] node1Id Id of the node at the other end of the pressure rate pump.
     * @param[in] pressure Pressure of the pump in Pa.
     * @return Id of the pressure pump.
     */
    int addPressurePump(int node0Id, int node1Id, double pressure);

    /**
     * @brief Specifies a node as sink.
     * @param[in] nodeId
     */
    void addSink(int nodeId);

    /**
     * @brief Adds or sets a node as the ground node, i.e., this node has a pressure value of 0 and acts as a reference node for all other nodes.
     * @param[in] nodeId
     */
    void addGround(int nodeId);

    /**
     * @brief Checks validity of a chip i.e. if network is one graph and all nodes and channels are connected to ground.
     * @return If chip is valid.
     */
    bool checkChipValidity();

    /**
     * @brief Add fluid to the simulator.
     * 
     * @param viscosity  Viscosity of the fluid in Pa s.
     * @param density Density of the fluid in kg/m^3.
     * @param concentration Concentration of the fluid in % (between 0.0 and 1.0).
     * @param molecularSize Molecular size in m^3.
     * @param diffusionCoefficient Diffusion coefficient of the fluid in m^2/s.
     * @param saturation Saturation value to translate the concentration in an actual concentration value in mol/m^3.
     * @return Id of the fluid.
     */
    int addFluid(double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient = 0.0, double saturation = 0.0);

    /**
     * @brief Add a mixture of fluids
     * 
     * @param[in] fluids a unordered map of fluid and concentration pairs.
     * @return Id of the fluid.
     */
    int addMixture(std::unordered_map<int, double> fluids);

    /**
     * @brief Specifies which fluid is the continuous phase.
     * @param[in] fluidId Id of the fluid that should be set as continuous phase.
     */
    void setContinuousPhase(int fluidId);

    /**
     * @brief Specifies the change of an input fluid from a certain pump starting at a certain time.
     * @param[in] fluidId Id of the fluid that should be injected.
     * @param[in] pumpId Id of the pump at which the fluid should be injected.
     * @param[in] injectionTime Time of the injection.
     */
    void setChangeInputFluid(int fluidId, int pumpId, double injectionTime);

    /**
     * @brief Specifies the change of an input mixture from a certain pump starting at a certain time.
     * 
     * @param mixtureId Id of the mixture that should be injected.
     * @param pumpId Id of the pump at which the fluid should be injected.
     * @param injectionTime Time of the injection.
     */
    void setChangeInputMixture(int mixtureId, int pumpId, double injectionTime);

    /**
     * @brief Define the maximal adaptive time step of the simulation.
    This time step is applied when a droplet changes channels in order to increase the simulation accuracy.
    A value of 0 disables this behavior (default is 0).
     * @param[in] timeStep in s.
     */
    void setMaximalAdaptiveTimeStep(double timeStep);

    /**
     * @brief Creates and adds a droplet to the simulation.
     * @param[in] fluidId Id of the fluid the droplet consists of.
     * @param[in] volume Volume of the droplet in m^3.
     * @param[in] injectionTime Simulation time at which the droplet should be injected.
     * @param[in] channelId Id of the channel in which the droplet is injected.
     * @param[in] relInjectionPosition Relative injection position (between 0.0 and 1.0) within the injection channel.
     * @return Id of the droplet.
     */
    int addDroplet(int fluidId, double volume, double injectionTime, int channelId, double relInjectionPosition);

    /**
     * @brief Set the duration which should be simulated. This is required for continuous fluid simulations.
     * 
     * @param duration Duration of the simulation in s.
     */
    void setSimulationDuration(double duration);

    /**
     * @brief Set the time step at which results should be stored. This is required for continuous fluid simulations. More timeSteps might happen in the background if necessary.
     * 
     * @param timeStep Time step in s.
     */
    void setSimulationResultTimeStep(double timeStep);

    /**
     * @brief Set the minimal internal simulation timestep for which an intermediate simulation state should be calculated. Required for continuous fluid simulations.
     * 
     * @param timeStep Minimal time interval for which the simulation should be calculated. More timeSteps might happen in the background if necessary.
     */
    void setSimulationCalculationTimeStep(double timeStep);

    /**
     * @brief Conduct the simulation.
     * @return The simulation result.
     */
    SimulationResult simulate();
};

}  // namespace droplet