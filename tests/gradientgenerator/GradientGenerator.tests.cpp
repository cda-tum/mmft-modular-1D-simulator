#include "Results.h"

#include <exception>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Simulator.h"
#include "gtest/gtest.h"

#include <execinfo.h>

/*
Continuous Fluid Simulation Tests based on
Gerold Fink et al. “Automatic Design of Microfluidic Gradient Generators”. In: IEEE Access 10 (2022), pp. 28155–28164. DOI: 10.1109/ACCESS.2022.3158327 (cit. on pp. 50, 51, 62).
*/

TEST(GradientGenerator, GradientGeneratorSmall) {  // Result 50%
    /* Settings Gradient Generator to generate this network:
    let w = 300e-6;
    let h = 100e-6;
    let radius = w;
    let wMeanderMax = 19 * w;
    let lConnection = (wMeanderMax + w) / 2;
    const lInlet = 10 * w;
    const lOutlet = 10 * w;
    let qInlet = 0.02 * w * h;
    let mu = 1e-3;
    let tMin = 1e-10;
    let concentrations = [0.5];
    let nOutlets = concentrations.length + 2;
    let nInlets = 2;
    */
    // std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    //std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = -6e-10;
    // create flow-rate pump from node -1 to node 0 with the given flow-rate
    auto pump0 = sim.addFlowRatePump(0, -1, flowRate);
    auto pump1 = sim.addFlowRatePump(1, -1, flowRate);

    //std::cout << "--- channels ---" << std::endl;
    // channels
    auto cWidth = 300e-6;
    auto cHeight = 100e-6;
    auto inletLength = 10 * cWidth;
    auto outletLength = 10 * cWidth;
    auto wMeanderMax = 19 * cWidth;
    auto connectionLength = (wMeanderMax + cWidth) / 2;

    // create channels
    // inlets
    auto c2 = sim.addChannel(0, 3, cHeight, cWidth, inletLength);
    auto c3 = sim.addChannel(1, 5, cHeight, cWidth, inletLength);
    // connections
    auto c4 = sim.addChannel(3, 2, cHeight, cWidth, connectionLength);
    auto c5 = sim.addChannel(3, 4, cHeight, cWidth, connectionLength);
    auto c6 = sim.addChannel(5, 4, cHeight, cWidth, connectionLength);
    auto c7 = sim.addChannel(5, 6, cHeight, cWidth, connectionLength);
    // meanders
    auto c8 = sim.addChannel(2, 7, cHeight, cWidth, 0.003084956172616956);
    auto c9 = sim.addChannel(4, 8, cHeight, cWidth, 0.004584956172616956);
    auto c10 = sim.addChannel(6, 9, cHeight, cWidth, 0.003084956172616958);
    // outlets
    auto c11 = sim.addChannel(7, -1, cHeight, cWidth, outletLength);
    auto c12 = sim.addChannel(8, -1, cHeight, cWidth, outletLength);
    auto c13 = sim.addChannel(9, -1, cHeight, cWidth, outletLength);

    //std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    //std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    //std::cout << "--- fluids ---" << std::endl;
    auto waterYellow = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    auto waterBlue = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    // auto water = sim.addFluid(8.65269e-4, 1.56e3, 1.0, 9e-10);  // 27 Grad
    // auto glucose = sim.addFluid(1.306489e-3, 1e3, 1.0, 9e-10);  // 27 Grad
    //std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(waterYellow);

    // continuous fluid injection
    //std::cout << "--- continuous injection ---" << std::endl;
    sim.setChangeInputFluid(waterBlue, pump1, 0.0);

    sim.setSimulationDuration(100.0);
    sim.setSimulationResultTimeStep(100.0);

    //std::cout << "--- validity check chip ---" << std::endl;

    // check if chip is valid
    sim.checkChipValidity();

    //std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network
    auto result = sim.simulate();

    // print the result
    //std::cout << "--- result ---" << std::endl;
    //std::cout << result.toJson(4) << std::endl;
    std::ofstream file("GradientGeneratorSmall.json");
    file << result.toJson(4);

    int lastStateId = result.states.back().id;
    auto fluidConcentrations11 = result.getAverageFluidConcentrationsInEdge(lastStateId, c11);
    ASSERT_NEAR(fluidConcentrations11.at(waterYellow), 1.0, std::numeric_limits<double>::epsilon());
    auto fluidConcentrations12 = result.getAverageFluidConcentrationsInEdge(lastStateId, c12);
    ASSERT_NEAR(fluidConcentrations12.at(waterYellow), 0.5, 5e-14);
    ASSERT_NEAR(fluidConcentrations12.at(waterBlue), 0.5, 5e-14);
    auto fluidConcentrations13 = result.getAverageFluidConcentrationsInEdge(lastStateId, c13);
    ASSERT_NEAR(fluidConcentrations13.at(waterBlue), 1.0, std::numeric_limits<double>::epsilon());
}

TEST(GradientGenerator, GradientGeneratorSmallDifferentPaper) {  // Result Paper 38.63%
    /* Settings Gradient Generator to generate this network: 
    let w = 300e-6; // mentioned in paper
    let h = 200e-6; // mentioned in paper
    let radius = w;
    let wMeanderMax = 21 * w;
    let lConnection = (wMeanderMax + w) / 2;
    const lInlet = 10 * w;
    const lOutlet = 10 * w;
    let qInlet = 5e-9 / 60; // mentioned in paper: 5uLsˆ-1
    let mu = 1e-3; 
    let tMin = 60;
    let concentrations = [0.3863]; // mentioned in paper
    let nOutlets = concentrations.length + 2; // due to concentrations
    let nInlets = 2; // due to concentrations
    */
    //std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    //std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = -5e-9 / 60;
    // create flow-rate pump from node -1 to node 0 with the given flow-rate
    auto pump0 = sim.addFlowRatePump(0, -1, flowRate);
    auto pump1 = sim.addFlowRatePump(1, -1, flowRate);

    //std::cout << "--- channels ---" << std::endl;
    // channels
    auto cWidth = 300e-6;
    auto cHeight = 200e-6;
    auto inletLength = 10 * cWidth;
    auto outletLength = 10 * cWidth;
    auto wMeanderMax = 21 * cWidth;
    auto connectionLength = (wMeanderMax + cWidth) / 2;

    // create channels
    // inlets
    auto c2 = sim.addChannel(0, 3, cHeight, cWidth, inletLength);
    auto c3 = sim.addChannel(1, 5, cHeight, cWidth, inletLength);
    // connections
    auto c4 = sim.addChannel(3, 2, cHeight, cWidth, connectionLength);
    auto c5 = sim.addChannel(3, 4, cHeight, cWidth, connectionLength);
    auto c6 = sim.addChannel(5, 4, cHeight, cWidth, connectionLength);
    auto c7 = sim.addChannel(5, 6, cHeight, cWidth, connectionLength);
    // meanders
    auto c8 = sim.addChannel(2, 7, cHeight, cWidth, 0.04163712059499694);
    auto c9 = sim.addChannel(4, 8, cHeight, cWidth, 0.055556211206648096);
    auto c10 = sim.addChannel(6, 9, cHeight, cWidth, 0.055783776606846545);
    // outlets
    auto c11 = sim.addChannel(7, -1, cHeight, cWidth, 0.006899999999999999);
    auto c12 = sim.addChannel(8, -1, cHeight, cWidth, 0.0009);
    auto c13 = sim.addChannel(9, -1, cHeight, cWidth, 0.006899999999999999);

    //std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    //std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    //std::cout << "--- fluids ---" << std::endl;
    auto waterYellow = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    auto waterBlue = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    //std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(waterYellow);

    // continuous fluid injection
    //std::cout << "--- continuous injection ---" << std::endl;
    sim.setChangeInputFluid(waterBlue, pump1, 0.0);

    sim.setSimulationDuration(85.0);
    sim.setSimulationResultTimeStep(85.0);

    //std::cout << "--- validity check chip ---" << std::endl;

    // check if chip is valid
    sim.checkChipValidity();

    //std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network

    auto result = sim.simulate();

    // print the result
    //std::cout << "--- result ---" << std::endl;
    // std::cout << result.toJson(4) << std::endl;
    std::ofstream file("GradientGeneratorSmallDifferent.json");
    file << result.toJson(4);

    int lastStateId = result.states.back().id;
    auto fluidConcentrations11 = result.getAverageFluidConcentrationsInEdge(lastStateId, c11);
    ASSERT_NEAR(fluidConcentrations11.at(waterYellow), 1.0, std::numeric_limits<double>::epsilon());
    auto fluidConcentrations12 = result.getAverageFluidConcentrationsInEdge(lastStateId, c12);
    ASSERT_NEAR(fluidConcentrations12.at(waterYellow), 0.3863, 5e-14);
    ASSERT_NEAR(fluidConcentrations12.at(waterBlue), 1.0 - 0.3863, 5e-14);
    auto fluidConcentrations13 = result.getAverageFluidConcentrationsInEdge(lastStateId, c13);
    ASSERT_NEAR(fluidConcentrations13.at(waterBlue), 1.0, std::numeric_limits<double>::epsilon());
}

TEST(GradientGenerator, GradientGeneratorUltraLargePaper) {  // Paper 100%/88.64%/40.12%/21.82%/0%
    /* Settings Gradient Generator to generate this network:
    let w = 300e-6;
    let h = 200e-6;
    let radius = w;
    let wMeanderMax = 21 * w;
    let lConnection = (wMeanderMax + w) / 2;
    const lInlet = 10 * w;
    const lOutlet = 10 * w;
    const straightOutlets = false;
    let qInlet = 5e-9 / 60;
    let mu = 1e-3;
    let tMin = 60;
    let concentrations = [0.8864, 0.4012, 0.2182];
    let nOutlets = concentrations.length + 2;
    let nInlets = 2;
    */
    //std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    //std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = -5e-9 / 60;
    // create flow-rate pump from node -1 to node 0 with the given flow-rate
    auto pump0 = sim.addFlowRatePump(0, -1, flowRate);
    auto pump1 = sim.addFlowRatePump(1, -1, flowRate);

    //std::cout << "--- channels ---" << std::endl;
    // channels
    auto cWidth = 300e-6;
    auto cHeight = 200e-6;
    auto inletLength = 10 * cWidth;
    auto outletLength = 10 * cWidth;
    auto wMeanderMax = 21 * cWidth;
    auto connectionLength = (wMeanderMax + cWidth) / 2;

    // create channels
    // inlets
    auto c2 = sim.addChannel(0, 3, cHeight, cWidth, inletLength);
    auto c3 = sim.addChannel(1, 5, cHeight, cWidth, inletLength);
    // connections
    auto c4 = sim.addChannel(3, 2, cHeight, cWidth, connectionLength);
    auto c5 = sim.addChannel(3, 4, cHeight, cWidth, connectionLength);
    auto c6 = sim.addChannel(5, 4, cHeight, cWidth, connectionLength);
    auto c7 = sim.addChannel(5, 6, cHeight, cWidth, connectionLength);
    // meanders
    auto c8 = sim.addChannel(2, 8, cHeight, cWidth, 0.0299671361633416);
    auto c9 = sim.addChannel(4, 10, cHeight, cWidth, 0.044431241343332664);
    auto c10 = sim.addChannel(6, 12, cHeight, cWidth, 0.03161589671051416);
    // connections
    auto c11 = sim.addChannel(8, 7, cHeight, cWidth, connectionLength);
    auto c12 = sim.addChannel(8, 9, cHeight, cWidth, connectionLength);
    auto c13 = sim.addChannel(10, 9, cHeight, cWidth, connectionLength);
    auto c14 = sim.addChannel(10, 11, cHeight, cWidth, connectionLength);
    auto c15 = sim.addChannel(12, 11, cHeight, cWidth, connectionLength);
    auto c16 = sim.addChannel(12, 13, cHeight, cWidth, connectionLength);
    // meanders
    auto c17 = sim.addChannel(7, 14, cHeight, cWidth, 0.03901274147675735);
    auto c18 = sim.addChannel(9, 15, cHeight, cWidth, 0.12054175031249897);
    auto c19 = sim.addChannel(11, 16, cHeight, cWidth, 0.047689442515671886);
    auto c20 = sim.addChannel(13, 17, cHeight, cWidth, 0.051441945635952405);
    // connections
    auto c21 = sim.addChannel(14, 18, cHeight, cWidth, connectionLength);
    auto c22 = sim.addChannel(14, 19, cHeight, cWidth, connectionLength);
    auto c23 = sim.addChannel(15, 19, cHeight, cWidth, connectionLength);
    auto c24 = sim.addChannel(15, 20, cHeight, cWidth, connectionLength);
    auto c25 = sim.addChannel(16, 20, cHeight, cWidth, connectionLength);
    auto c26 = sim.addChannel(16, 21, cHeight, cWidth, connectionLength);
    auto c27 = sim.addChannel(17, 21, cHeight, cWidth, connectionLength);
    auto c28 = sim.addChannel(17, 22, cHeight, cWidth, connectionLength);
    // meanders
    auto c29 = sim.addChannel(18, 23, cHeight, cWidth, 0.026532242716654282);
    auto c30 = sim.addChannel(19, 24, cHeight, cWidth, 0.03333335815696055);
    auto c31 = sim.addChannel(20, 25, cHeight, cWidth, 0.03948202927907311);
    auto c32 = sim.addChannel(21, 26, cHeight, cWidth, 0.03355323373613861);
    auto c33 = sim.addChannel(22, 27, cHeight, cWidth, 0.02498585157999049);
    // outlets
    auto c34 = sim.addChannel(23, -1, cHeight, cWidth, 0.013499999999999998);
    auto c35 = sim.addChannel(24, -1, cHeight, cWidth, 0.007499999999999999);
    auto c36 = sim.addChannel(25, -1, cHeight, cWidth, 0.0014999999999999998);
    auto c37 = sim.addChannel(26, -1, cHeight, cWidth, 0.007499999999999999);
    auto c38 = sim.addChannel(27, -1, cHeight, cWidth, 0.013499999999999998);

    //std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    //std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    //std::cout << "--- fluids ---" << std::endl;
    // auto water = sim.addFluid(8.65269e-4, 1.56e3, 9e-10, 0.0);  // 27 Grad
    // //auto glucose = sim.addFluid(1.306489e-3, 1e3, 9e-10, 1.0);  // 27 Grad
    // auto glucose = sim.addFluid(8.65269e-4, 1.56e3, 9e-10, 1.0);  // 27 Grad
    auto waterYellow = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    auto waterBlue = sim.addFluid(1e-3, 1.56e3, 1.0, 9e-10);
    //std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(waterYellow);

    // continuous fluid injection
    //std::cout << "--- continuous injection ---" << std::endl;
    sim.setChangeInputFluid(waterBlue, pump1, 0.0);

    // sim.setSimulationDuration(560.0);
    // sim.setSimulationResultTimeStep(560.0);
    sim.setSimulationDuration(1000.0);
    sim.setSimulationResultTimeStep(1000.0);

    //std::cout << "--- validity check chip ---" << std::endl;

    // check if chip is valid
    sim.checkChipValidity();

    //std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network

    auto result = sim.simulate();

    // print the result
    //std::cout << "--- result ---" << std::endl;
    // std::cout << result.toJson(4) << std::endl;
    std::ofstream file("GradientGeneratorUltraLarge.json");
    file << result.toJson(4);

    // Paper 100%/88.64%/40.12%/21.82%/0%
    int lastStateId = result.states.back().id;
    auto fluidConcentrations34 = result.getAverageFluidConcentrationsInEdge(lastStateId, c34);
    ASSERT_NEAR(fluidConcentrations34.at(waterYellow), 1.0, std::numeric_limits<double>::epsilon());
    auto fluidConcentrations35 = result.getAverageFluidConcentrationsInEdge(lastStateId, c35);
    ASSERT_NEAR(fluidConcentrations35.at(waterYellow), 0.8864, 5e-14);
    ASSERT_NEAR(fluidConcentrations35.at(waterBlue), 1 - 0.8864, 5e-14);
    auto fluidConcentrations36 = result.getAverageFluidConcentrationsInEdge(lastStateId, c36);
    ASSERT_NEAR(fluidConcentrations36.at(waterYellow), 0.4012, 5e-14);
    ASSERT_NEAR(fluidConcentrations36.at(waterBlue), 1.0 - 0.4012, 5e-14);
    auto fluidConcentrations37 = result.getAverageFluidConcentrationsInEdge(lastStateId, c37);
    ASSERT_NEAR(fluidConcentrations37.at(waterYellow), 0.2182, 5e-14);
    ASSERT_NEAR(fluidConcentrations37.at(waterBlue), 1.0 - 0.2182, 5e-14);
    auto fluidConcentrations38 = result.getAverageFluidConcentrationsInEdge(lastStateId, c38);
    ASSERT_NEAR(fluidConcentrations38.at(waterBlue), 1.0, std::numeric_limits<double>::epsilon());
}
