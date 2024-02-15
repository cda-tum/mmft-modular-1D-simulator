/**
 * @file Simulation.h
 */

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include "ChangeInputMixture.h"
#include "Chip.h"
#include "Droplet.h"
#include "DropletBoundary.h"
#include "Event.h"
#include "Fluid.h"
#include "IMembraneResistanceModel.h"
#include "IResistanceModel.h"
#include "Injection.h"
#include "MembraneResistanceModels.h"
#include "Mixture.h"
#include "Results.h"

namespace sim {

/**
 * @brief Enum to define the available resistance models.
 */
enum class ResistanceModel {
    ONE_D_MODEL,             ///< 1D Resistance Model.
    TEST_MODEL,              ///< Test Resistance Model.
    ONE_D_CONTINUOUS_MODEL,  ///< 1D Resistance Model.
};

/**
 * @brief Enum to define the available membrane resistance models.
 */
enum class MembraneResistanceModel { MODEL_0, MODEL_1, MODEL_2, MODEL_3, MODEL_4, MODEL_5, MODEL_6, MODEL_7, MODEL_8, MODEL_9 };

/**
 * @brief Class that conducts the simulation and owns all parameters necessary for it.
 */
class Simulation {
  private:
    arch::Chip* chip;                                                                  ///< Chip for which the simulation should be conducted.
    std::unordered_map<int, std::unique_ptr<Fluid>> fluids;                            ///< Fluids specified for the simulation.
    std::unordered_map<int, std::unique_ptr<Droplet>> droplets;                        ///< Droplets which are simulated.
    std::unordered_map<int, std::unique_ptr<Injection>> injections;                    ///< Injections of droplets that should take place during the simulation.
    std::unordered_map<int, std::unique_ptr<ChangeInputMixture>> changeInputMixtures;  ///< Injections of fluids that should take place during the simulation.
    std::unordered_map<int, std::deque<std::pair<int, double>>> mixturesInEdge;        ///< Which mixture currently flows in which edge <EdgeID, <MixtureID, currPos>>>
    std::vector<Mixture> mixtures;
    // std::unordered_map<int, std::unique_ptr<Mixture>> mixtures;                        ///< Fluid mixtures, contains fluidIds and concentrations of which a mixture is made of
    std::unique_ptr<IResistanceModel> resistanceModel;                                       ///< The resistance model used for the simulation.
    std::unique_ptr<IMembraneResistanceModel> membraneResistanceModel;                       ///< The membrane resistance model used for the simulation.
    ResistanceModel resistanceModelName = ResistanceModel::ONE_D_MODEL;                      ///< Which resistance model should be used for the calculations of this simulation.
    MembraneResistanceModel membraneResistanceModelName = MembraneResistanceModel::MODEL_9;  ///< Which membrane resistance model should be used for the calculations of this simulation.
    int continuousPhaseMixtureId = -1;                                                       ///< MixtureId of the continuous phase.
    double const slipFactor = 1.28;                                                          ///< Slip factor of droplets.
    double maximalAdaptiveTimeStep = 0;                                                      ///< Maximal adaptive time step that is applied when droplets change the channel.
    double currTime = 0;                                                                     ///< The time elapsed since the start of the simulation in s.
    double iState = 0;                                                                       ///< The current index of the next state.
    bool dropletsAtBifurcation = false;                                                      ///< If one or more droplets are currently at a bifurcation. Triggers the usage of the maximal adaptive time step.
    bool enableMerging = true;                                                               ///< If the droplet merging simulation should be performed.
    bool simulationInProgress;                                                               ///< Information if the simulation is currently in progress.
    double simulationDuration = std::numeric_limits<double>::max();                          ///< Simulation duration in s. Required for continuous fluid simulations.
    double simulationResultTimeStep = 0.2;                                                   ///< Simulation time steps for the output of the simulation in s. Required for continuous fluid simulations.
    double minSimulationCalcTimeStep = std::numeric_limits<double>::max();                   ///< Minimal simulation time steps for the internal simulation calculation.
    double internalSimulationTimeStep;                                                       ///< Internal simulation time steps for which the new state of the simulation is calculated (if no other events take place). Required for continuous fluid simulations. Might be higher than simulationTimeStep and is based on the shortest throughput time of an edge in the network.

    /**
     * @brief Initializes the resistance model and the channel resistances of the empty channels.
     */
    void initialize();

    /**
     * @brief Get the average viscosity of mixtures in a channel 
     * 
     * @param channelId Id to channel for which average viscosity should be calculated
     * @return double Average viscosity of mixtures in channel
     */
    double getAverageViscosityInChannel(int channelId);

    /**
     * @brief Update the channel resistances based on the current average viscosity of fluids within a channel. (Only has an effect if the 1D continuous resistance model is used, as otherwise the continuous phase viscosity is considered which does not change.)
     */
    void updateChannelResistances();

    /**
     * @brief Update the droplet resistances of the channels based on the current positions of the droplets.
     */
    void updateDropletResistances();

    /**
     * @brief Compute all possible next events.
     */
    std::vector<std::unique_ptr<Event>> computeEvents(double simulationResultTimeCounter);

    /**
     * @brief Moves all droplets according to the given time step.
     * @param[in] timeStep Duration to which the droplet movement should be forwarded.
     */
    void moveDroplets(double timeStep);

    /**
     * @brief Get the outflow volume of a pump at a node.
     * 
     * @param nodeId Id of the node.
     * @param flowRate Flow-rate of the pump.
     * @param timeStep Time step for which the outflow volume should be calculated.
     * @return Outflow volume for a specific duration at a specific node.
     */
    double getPumpOutflowVolume(int nodeId, double flowRate, double timeStep);

    /**
     * @brief Calculate and set new state of the continuous fluid simulation. Move mixture positions and create new mixtures if necessary.
     * 
     * @param timeStep Time step in s for which the new mixtures state should be calculated.
     */
    void calculateNewMixtures(double timeStep);

    /**
     * @brief Store simulation parameters to the result.
     * @param[in, out] result Reference to the simulation result in which all current parameters of the simulation should be stored.
     */
    void storeSimulationParameters(droplet::SimulationResult& result);

    /**
     * @brief Store all simulation results to the result.
     * @param[in, out] result Reference to the simulation result in which all current parameters of the simulation should be stored.
     */
    void storeSimulationResults(droplet::SimulationResult& result);

  public:
    /**
     * @brief Creates simulation.
     */
    Simulation();

    /**
     * @brief Set the chip for which the simulation should be conducted.
     * @param[in] chip Chip on which the simulation will be conducted.
     */
    void setChip(arch::Chip* chip);

    /**
     * @brief Get the chip.
     * @return Chip or nullptr if no chip is specified.
     */
    arch::Chip* getChip();

    /**
     * @brief Define which fluid should act as continuous phase, i.e., as carrier fluid for the droplets.
     * @param[in] fluidId Unique identifier of the fluid the continuous phase consists of.
     */
    void setContinuousPhase(int fluidId);

    /**
     * @brief Get the Mixture that consists of 100% of the fiven fluid object.
     * 
     * @param fluidId Unique identifier of the fluid.
     */
    int getMixtureOfFluid(int fluidId);

    /**
     * @brief Fills the chip with the continuous phase mixture (in the beginning of the simulation).
     * 
     */
    void fillChipWithContinuousPhaseMixture();

    /**
     * @brief Set a pointer to the Mixture Injection object.
     * 
     * @param fluidId Unique identifier of the fluid to be injected.
     * @param pumpId Unique identifier of the pump at which the fluid should be injected.
     * @param injectionTime Time at which the fluid injection should start.
     * @return Pointer to created injection.
     */
    ChangeInputMixture* setChangeInputFluid(int fluidId, int pumpId, double injectionTime);

    /**
     * @brief Set a pointer to the Mixture Injection object.
     * 
     * @param mixtureId Unique identifier of the mixture to be injected.
     * @param pumpId Unique identifier of the pump at which the fluid should be injected.
     * @param injectionTime Time at which the fluid injection should start.
     * @return Pointer to created injection.
     */
    ChangeInputMixture* setChangeInputMixture(int mixtureId, int pumpId, double injectionTime);

    /**
     * @brief Get a pointer to the Fluid Injection object.
     * 
     * @param injectionId Id of the continous injection.
     * @return Pointer to the continuous injection.
     */
    ChangeInputMixture* getChangeInputFluid(int injectionId);

    /**
     * @brief Define which resistance model should be used for the channel and droplet resistance calculations.
     * @param[in] modelName Name of the resistance model to be used.
     */
    void setResistanceModel(ResistanceModel modelName);

    /**
     * @brief Set the membrane resistance model that should be used for the membrane resistance calculation during the simulation.
     * 
     * @param modelName Name of the membrane resistance model.
     */
    void setMembraneResistanceModel(MembraneResistanceModel modelName);

    /**
     * @brief Define the maximal adaptive time step of the simulation.
     This time step is applied when a droplet changes channels in order to increase the simulation accuracy.
     A value of 0 disables this behavior (default is 0).
     * @param[in] timeStep Maximal time step that the simulation can do when changing a channel.
     */
    void setMaximalAdaptiveTimeStep(double timeStep);

    /**
     * @brief Add new fluid to the simulation.
     * 
     * @param viscosity Dynamic viscosity of the fluid in Pas.
     * @param density Density of the fluid in kg/m^3.
     * @param concentration Concentration of the fluid in percent (between 0.0 and 1.0).
     * @param molecularSize Molecular size in m^3.
     * @param diffusionCoefficient Diffusion coefficient of the fluid in m^2/s.
     * @param saturation Saturation value to translate the concentration in an actual concentration value [mol/m^3].
     * @return Pointer to the fluid.
     */
    Fluid* addFluid(double viscosity, double density, double concentration, double molecularSize, double diffusionCoefficient, double saturation);

    /**
     * @brief Get fluid.
     * @param[in] fluidId Id of the fluid
     * @return Pointer to fluid with the corresponding id
     */
    Fluid* getFluid(int fluidId);

    /**
     * @brief Get mixture.
     * 
     * @param mixtureId Id of the mixture
     * @return Pointer to mixture with the correspondig id
     */
    Mixture* getMixture(int mixtureId);

    /**
     * @brief Add mixture
     *
     * @param[in] fluids List of fluids and concentration pairs. 
     */
    int addMixture(std::unordered_map<int, double>&& fluidConcentrations, int id = -1);

    /**
     * @brief Create droplet.
     * @param[in] fluidId Unique identifier of the fluid the droplet consists of.
     * @param[in] volume Volume of the fluid in m^3.
     * @return Pointer to created droplet.
     */
    Droplet* addDroplet(int fluidId, double volume);

    /**
     * @brief Get droplet
     * @param dropletId Id of the droplet
     * @return Pointer to droplet with the corresponding id
     */
    Droplet* getDroplet(int dropletId);

    /**
     * @brief Gets droplet that is present at the corresponding node (i.e., the droplet spans over this node).
     * @param nodeId The id of the node
     * @return Pointer to droplet or nullptr if no droplet was found
     */
    Droplet* getDropletAtNode(int nodeId);

    /**
     * @brief Create injection.
     * @param[in] dropletId Id of the droplet that should be injected.
     * @param[in] injectionTime Time at which the droplet should be injected in s.
     * @param[in] channelId Id of the channel, where droplet should be injected.
     * @param[in] injectionPosition Position inside the channel at which the droplet should be injected (relative to the channel length between 0.0 and 1.0).
     * @return Pointer to created injection.
     */
    Injection* addInjection(int dropletId, double injectionTime, int channelId, double injectionPosition);

    /**
     * @brief Get injection
     * @param injectionId The id of the injection
     * @return Pointer to injection with the corresponding id.
     */
    Injection* getInjection(int injectionId);

    /**
     * @brief Creates a new fluid out of two existing fluids.
     * @param fluid0Id Id of the first fluid.
     * @param volume0 The volume of the first fluid.
     * @param fluid1Id Id of the second fluid.
     * @param volume1 The volume of the second fluid.
     * @return Pointer to new fluid.
     */
    Fluid* mixFluids(int fluid0Id, double volume0, int fluid1Id, double volume1);

    /**
     * @brief Creates a new droplet from two existing droplets.
     Please note that this only creates a new droplet inside the simulation, but the actual boundaries have to be set separately, which is usually done inside the corresponding merge events.
     * @param droplet0Id Id of the first droplet.
     * @param droplet1Id Id of the second droplet.
     * @return Pointer to new droplet.
     */
    Droplet* mergeDroplets(int droplet0Id, int droplet1Id);

    /**
     * @brief Set the total duration of the simulation.
     * 
     * @param duration Total duration in [s].
     */
    void setSimulationDuration(double duration);

    /**
     * @brief Set the timestep for the result output. 
     * 
     * @param timeStep Time interval in [s].
     */
    void setSimulationResultTimeStep(double timeStep);

    /**
     * @brief Set the internal calulcation time step.
     * 
     * @param timeStep Time interval in [s].
     */
    void setSimulationCalculationTimeStep(double timeStep);

    /**
     * @brief Conduct the simulation.
     * @return The result of the simulation containing all intermediate simulation steps and calculated parameters.
     */
    droplet::SimulationResult simulate();

    /**
     * @brief Stop the simulation.
     * 
     */
    void stopSimulation();
};

}  // namespace sim