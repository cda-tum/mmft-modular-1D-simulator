#include <fstream>
#include <iostream>

#include "droplet-simulator/Results.h"
#include "droplet-simulator/Simulator.h"

int main(int argc, char const* argv[]) {
    std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = 3e-3;
    // create flow-rate pump from node -1 to node 0 with the given flow-rate
    auto pump = sim.addFlowRatePump(-1, 0, flowRate);

    std::cout << "--- channels ---" << std::endl;
    // channels
    auto cWidth = 100e-6;
    auto cHeight = 30e-6;
    auto cLength = 1000e-6;

    // create channel from node 0 to node 1 with the given height, width, length
    auto c1 = sim.addChannel(0, 1, cHeight, cWidth, cLength);
    // create channel from node 1 to node 2 with the given height, width, length
    auto c2 = sim.addChannel(1, 2, cHeight, cWidth, cLength);
    // create channel from node 2 to node 3 with the given height, width, length
    auto c3 = sim.addChannel(2, 3, cHeight, cWidth, cLength);
    // create channel from node 2 to node 4 with the given height, width, length
    auto c4 = sim.addChannel(3, 4, cHeight, cWidth, cLength);
    // create channel from node 3 to node 4 with the given height, width, length
    auto c5 = sim.addChannel(4, 5, cHeight, cWidth, cLength);
    // create channel from node 4 to node -1 with the given height, width, length
    auto c6 = sim.addChannel(5, -1, cHeight, cWidth, cLength);
    std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    std::cout << "--- fluids ---" << std::endl;
    // add fluid with 1e-3 viscosity and 1e3 density
    auto fluid0 = sim.addFluid(1e-3, 1e3, 0.0, 9e-10);  // wasser mit 0% glucose
    // add fluid with 3e-3 viscosity and 1e3 density
    auto fluid1 = sim.addFluid(3e-3, 1e3, 1.0, 9e-10);  // 100% glucose
    auto fluid2 = sim.addFluid(3e-3, 1e3, 0.5, 9e-10);  // 50% glucose
    std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(fluid0);

    // droplet
    // std::cout << "--- droplet ---" << std::endl;
    // auto dropletVolume = 1.5 * cWidth * cWidth * cHeight;
    // // create a droplet of fluid1, with a given droplet volume, injected at injectionTime 0.0 in channel with the channelId c1 at relative injection position 0.5
    // auto droplet0 = sim.addDroplet(fluid1, dropletVolume, 0.0, c1, 0.5);

    // continuous fluid injection
    std::cout << "--- continuous injection ---" << std::endl;
    sim.setChangeInputFluid(fluid1, pump, 0.0);
    sim.setChangeInputFluid(fluid2, pump, 1.0);

    // membrane and organ
    std::cout << "--- membrane ---" << std::endl;
    auto membrane0 = sim.addMembraneToChannel(3, cHeight, cWidth, 1e-6, 0.5);
    auto organ0 = sim.addOrganToMembrane(membrane0, cHeight, cWidth);

    sim.setSimulationDuration(4.0);
    sim.setSimulationResultTimeStep(0.2);

    std::cout << "--- validity check chip ---" << std::endl;
    // check if chip is valid
    sim.checkChipValidity();

    std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network
    auto result = sim.simulate();

    // print the result
    std::cout << "--- result ---" << std::endl;
    std::cout << result.toJson(4) << std::endl;

    std::ofstream file("result.json");
    file << result.toJson(4);

    return 0;
}
