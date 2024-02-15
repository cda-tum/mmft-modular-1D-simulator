#include "Results.h"

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Simulator.h"
#include "gtest/gtest.h"

/*
Membrane and Organ Tank Testcase based on 

Alan Chramiec et al. “Integrated human organ-on-a-chip model for predictive studies of anti-tumor drug efficacy and cardiac safety”. In: Lab Chip 20 (23 2020), pp. 4357–4372. DOI: 10.1039/D0LC00424C. URL: http://dx.doi.org/10.1039/ D0LC00424C (cit. on pp. 32, 59, 60, 61).

executed with MembraneResistanceModel9
*/
TEST(Organ, TwoOrgan) {
    //std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    //std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = 5.5e-8;  //3.3mL/min
    // create flow-rate pump from node -1 to node 0 with the given flow-rate
    auto pump0 = sim.addFlowRatePump(0, -1, flowRate);

    //std::cout << "--- channels ---" << std::endl;
    // channels
    auto cWidth = 5e-3;
    auto cHeight = 0.3e-3;          // mentioned in paper
    auto cConnectionLength = 5e-3;  // mentioned in paper

    auto cMembraneLength = 8e-3;   // mentioned in paper
    auto mHeight = 55e-6;          // membrane as stated in personal communication with the authors and according to specification on website
    auto mWidth = 4e-3;            // mentioned in paper
    auto pillarArea = 6.28318e-6;  // area covered by pillars
    auto poreRadius = 20e-6 / 2;   // pore size 20e-6 mentioned in paper, millipore size specifies size as diameter
    auto porosity = 0.14;

    auto oHeight = 13e-3;
    auto oWidth = 1.5e-6 / (oHeight * cMembraneLength);  //1.5mL organ volume mentioned in paper

    // create channels
    auto c1 = sim.addChannel(0, 1, cHeight, cWidth, cMembraneLength);  //bone tumor tissue
    auto c2 = sim.addChannel(1, 2, cHeight, cWidth, cConnectionLength);
    auto c3 = sim.addChannel(2, -1, cHeight, cWidth, cMembraneLength);  //cardiac tissue
    // create membrane and organ
    auto m4 = sim.addMembraneToChannel(c1, mHeight, mWidth, poreRadius, porosity);
    auto o5 = sim.addOrganToMembrane(m4, oHeight, oWidth);

    auto m6 = sim.addMembraneToChannel(c3, mHeight, mWidth - pillarArea / cMembraneLength, poreRadius, porosity);
    auto o7 = sim.addOrganToMembrane(m6, oHeight, oWidth);

    //std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    //std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    //std::cout << "--- fluids ---" << std::endl;
    auto continuousPhaseFluid = sim.addFluid(0.7e-3, 0.993e3, 1.0, 9e-10, 4.4e-10, 0.0);
    auto injectionFluid = sim.addFluid(0.7e-3, 0.993e3, 1.0, 9e-10, 4.4e-10, 3.894e-3);  // Linsitinib

    //std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(continuousPhaseFluid);

    // continuous fluid injection
    //std::cout << "--- continuous injection ---" << std::endl;
    sim.setChangeInputFluid(injectionFluid, pump0, 0.0);

    sim.setSimulationDuration(86400.0);  //24H
    // sim.setSimulationTimeStep(7200.0);   //2H
    sim.setSimulationResultTimeStep(1800.0);  //0.5H
    //sim.setSimulationCalculationTimeStep(0.01);

    //std::cout << "--- validity check chip ---" << std::endl;

    // check if chip is valid
    sim.checkChipValidity();

    //std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network
    auto result = sim.simulate();

    // print the result
    //std::cout << "--- result ---" << std::endl;
    //std::cout << result.toJson(4) << std::endl;
    std::ofstream file("TwoOrgan.json");
    file << result.toJson(4);

    auto fluidConcentrations3 = result.getAverageFluidConcentrationsInEdge(3 / 0.5, o5);
    ASSERT_NEAR(fluidConcentrations3.at(injectionFluid), 0.7, 0.15);
    auto fluidConcentrations6 = result.getAverageFluidConcentrationsInEdge(6 / 0.5, o5);
    ASSERT_NEAR(fluidConcentrations6.at(injectionFluid), 0.9, 0.1);
    auto fluidConcentrations12 = result.getAverageFluidConcentrationsInEdge(12 / 0.5, o5);
    ASSERT_NEAR(fluidConcentrations12.at(injectionFluid), 0.96, 0.05);
    auto fluidConcentrations24 = result.getAverageFluidConcentrationsInEdge(24 / 0.5, o5);
    ASSERT_NEAR(fluidConcentrations24.at(injectionFluid), 1.0, 0.05);

    fluidConcentrations3 = result.getAverageFluidConcentrationsInEdge(3 / 0.5, o7);
    ASSERT_NEAR(fluidConcentrations3.at(injectionFluid), 0.6, 0.1);
    fluidConcentrations6 = result.getAverageFluidConcentrationsInEdge(6 / 0.5, o7);
    ASSERT_NEAR(fluidConcentrations6.at(injectionFluid), 0.8, 0.1);
    fluidConcentrations12 = result.getAverageFluidConcentrationsInEdge(12 / 0.5, o7);
    ASSERT_NEAR(fluidConcentrations12.at(injectionFluid), 0.92, 0.05);
    fluidConcentrations24 = result.getAverageFluidConcentrationsInEdge(24 / 0.5, o7);
    ASSERT_NEAR(fluidConcentrations24.at(injectionFluid), 1.0, 0.05);

    // auto tumorConcentration = 0.0;
    // auto cardiacConcentration = 0.0;
    // auto channelConcentration = 0.0;
    // for (auto& state : result.states) {
    //     std::cout << "state " << state.id << " hour " << state.id * 0.5 << std::endl;
    //     auto averageFluidConcentrationsInEdge = result.getAverageFluidConcentrationsInEdge(state.id, o5);
    //     auto tumorConcentrationIncrease = 0.0;
    //     if (averageFluidConcentrationsInEdge.count(injectionFluid) > 0) {
    //         tumorConcentrationIncrease = averageFluidConcentrationsInEdge.at(injectionFluid) - tumorConcentration;
    //         tumorConcentration = averageFluidConcentrationsInEdge.at(injectionFluid);
    //     }
    //     std::cout << std::setprecision(15) << "tumor ";
    //     std::cout << " Concentration [mM] " << tumorConcentration;
    //     std::cout << " Increase " << tumorConcentrationIncrease << std::endl;

    //     averageFluidConcentrationsInEdge = result.getAverageFluidConcentrationsInEdge(state.id, o7);
    //     auto cardiacConcentrationIncrease = 0.0;
    //     if (averageFluidConcentrationsInEdge.count(injectionFluid) > 0) {
    //         cardiacConcentrationIncrease = averageFluidConcentrationsInEdge.at(injectionFluid) - cardiacConcentration;
    //         cardiacConcentration = averageFluidConcentrationsInEdge.at(injectionFluid);
    //     }
    //     std::cout << "cardiac ";
    //     std::cout << " Concentration [mM] " << cardiacConcentration;
    //     std::cout << " Increase " << cardiacConcentrationIncrease << std::endl;

    //     std::cout << "tumor - cardiac " << tumorConcentration - cardiacConcentration << std::endl;

    //     averageFluidConcentrationsInEdge = result.getAverageFluidConcentrationsInEdge(state.id, c2);
    //     auto channelConcentrationIncrease = 0.0;
    //     if (averageFluidConcentrationsInEdge.count(injectionFluid) > 0) {
    //         cardiacConcentrationIncrease = averageFluidConcentrationsInEdge.at(injectionFluid) - channelConcentration;
    //         channelConcentration = averageFluidConcentrationsInEdge.at(injectionFluid);
    //     }
    //     std::cout << "channel2 ";
    //     std::cout << " Concentration [mM] " << channelConcentration;
    //     std::cout << " Increase " << channelConcentrationIncrease << std::endl;

    //     std::cout << std::endl;
    // }
}
