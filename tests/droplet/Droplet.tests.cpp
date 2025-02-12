#include "Results.h"

#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Simulator.h"
#include "gtest/gtest.h"

/*
Droplet Simulation Test based on 
Gerold Fink et al. “Automatic Design of Droplet-Based Microfluidic Ring Networks”. In: IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems 40.2 (2021), pp. 339–349. DOI: 10.1109/TCAD.2020.2997000 (cit. on pp. 27, 54, 55, 56, 57).
*/

TEST(Droplet, RingNetworkE1) {
    //--- API ---
    droplet::Simulator sim;

    // flowRate pump
    auto flowRate = 3e-11;
    auto pump = sim.addFlowRatePump(-1, 0, flowRate);

    // droplet length
    auto lDroplet = 150e-6;

    // droplet distances
    double d_m1 = 881e-6;
    double d_m2 = 41230e-6;
    double d_m3 = 41867e-6;
    double d_m4 = 41956e-6;
    double d_m5 = 42047e-6;
    double d_m6 = 42127e-6;

    // channels
    auto cWidth = 100e-6;
    auto cHeight = 33e-6;
    auto bWidth = 200e-6;

    auto iLength = 2000e-6;

    // Node 1
    auto lIn1N1 = 469e-6;
    auto lIn2N1 = 531e-6;
    auto lByN1 = 1200e-6;
    auto lOut1N1 = 14531e-6;
    auto lOut2N1 = 24469e-6;
    // Node 2
    auto lIn1N2 = 743e-6;
    auto lIn2N2 = 757e-6;
    auto lByN2 = 1700e-6;
    auto lOut1N2 = 14257e-6;
    auto lOut2N2 = 24243e-6;
    // Node 3
    auto lIn1N3 = 1021e-6;
    auto lIn2N3 = 979e-6;
    auto lByN3 = 2200e-6;
    auto lOut1N3 = 13979e-6;
    auto lOut2N3 = 24021e-6;
    // Node 4
    auto lIn1N4 = 1300e-6;
    auto lIn2N4 = 1200e-6;
    auto lByN4 = 2700e-6;
    auto lOut1N4 = 13700e-6;
    auto lOut2N4 = 23800e-6;
    // Node 5
    auto lIn1N5 = 1580e-6;
    auto lIn2N5 = 1420e-6;
    auto lByN5 = 3200e-6;
    auto lOut1N5 = 13420e-6;
    auto lOut2N5 = 23580e-6;
    // Node 6
    auto lIn1N6 = 1862e-6;
    auto lIn2N6 = 1638e-6;
    auto lByN6 = 3700e-6;
    auto lOut1N6 = 13138e-6;
    auto lOut2N6 = 23362e-6;

    auto cLength = 500e-6;

    auto dLength1 = 881e-6;
    auto dLength2 = 40678e-6;
    auto dLength3 = 42185e-6;
    auto dLength4 = 41475e-6;

    auto c1 = sim.addChannel(0, 1, cHeight, cWidth, iLength);
    //Node 1
    auto c2_n1 = sim.addChannel(1, 2, cHeight, cWidth, lIn1N1);
    auto c3_n1 = sim.addChannel(1, 3, cHeight, cWidth, lIn2N1);
    auto c4_n1 = sim.addChannel(2, 4, cHeight, cWidth, lOut1N1);
    auto c5_n1 = sim.addChannel(3, 4, cHeight, cWidth, lOut2N1);
    auto c6_n1 = sim.addChannel(4, 21, cHeight, cWidth, cLength);
    auto c7_n1 = sim.addBypassChannel(2, 3, cHeight, bWidth, lByN1);
    //Node 2
    auto c2_n2 = sim.addChannel(21, 22, cHeight, cWidth, lIn1N2);
    auto c3_n2 = sim.addChannel(21, 23, cHeight, cWidth, lIn2N2);
    auto c4_n2 = sim.addChannel(22, 24, cHeight, cWidth, lOut1N2);
    auto c5_n2 = sim.addChannel(23, 24, cHeight, cWidth, lOut2N2);
    auto c6_n2 = sim.addChannel(24, 31, cHeight, cWidth, cLength);
    auto c7_n2 = sim.addBypassChannel(22, 23, cHeight, bWidth, lByN2);
    //Node 3
    auto c2_n3 = sim.addChannel(31, 32, cHeight, cWidth, lIn1N3);
    auto c3_n3 = sim.addChannel(31, 33, cHeight, cWidth, lIn2N3);
    auto c4_n3 = sim.addChannel(32, 34, cHeight, cWidth, lOut1N3);
    auto c5_n3 = sim.addChannel(33, 34, cHeight, cWidth, lOut2N3);
    auto c6_n3 = sim.addChannel(34, 41, cHeight, cWidth, cLength);
    auto c7_n3 = sim.addBypassChannel(32, 33, cHeight, bWidth, lByN3);
    //Node 4
    auto c2_n4 = sim.addChannel(41, 42, cHeight, cWidth, lIn1N4);
    auto c3_n4 = sim.addChannel(41, 43, cHeight, cWidth, lIn2N4);
    auto c4_n4 = sim.addChannel(42, 44, cHeight, cWidth, lOut1N4);
    auto c5_n4 = sim.addChannel(43, 44, cHeight, cWidth, lOut2N4);
    auto c6_n4 = sim.addChannel(44, 51, cHeight, cWidth, cLength);
    auto c7_n4 = sim.addBypassChannel(42, 43, cHeight, bWidth, lByN4);
    //Node 5
    auto c2_n5 = sim.addChannel(51, 52, cHeight, cWidth, lIn1N5);
    auto c3_n5 = sim.addChannel(51, 53, cHeight, cWidth, lIn2N5);
    auto c4_n5 = sim.addChannel(52, 54, cHeight, cWidth, lOut1N5);
    auto c5_n5 = sim.addChannel(53, 54, cHeight, cWidth, lOut2N5);
    auto c6_n5 = sim.addChannel(54, 61, cHeight, cWidth, cLength);
    auto c7_n5 = sim.addBypassChannel(52, 53, cHeight, bWidth, lByN5);
    //Node 6
    auto c2_n6 = sim.addChannel(61, 62, cHeight, cWidth, lIn1N6);
    auto c3_n6 = sim.addChannel(61, 63, cHeight, cWidth, lIn2N6);
    auto c4_n6 = sim.addChannel(62, 64, cHeight, cWidth, lOut1N6);
    auto c5_n6 = sim.addChannel(63, 64, cHeight, cWidth, lOut2N6);
    auto c6_n6 = sim.addChannel(64, -1, cHeight, cWidth, cLength);
    auto c7_n6 = sim.addBypassChannel(62, 63, cHeight, bWidth, lByN6);

    //--- sink ---
    sim.addSink(-1);
    //--- ground ---
    sim.addGround(-1);

    // fluids
    auto fluid0 = sim.addFluid(1e-3, 1e3, 0.0, 9e-10);
    auto fluid1 = sim.addFluid(4e-3, 1e3, 1.0, 9e-10);
    //--- continuousPhase ---
    sim.setContinuousPhase(fluid0);

    // droplets
    auto dropletVolume = lDroplet * cWidth * cHeight;
    auto droplet0 = sim.addDroplet(fluid1, dropletVolume, 0.0, c1, 0.5);
    auto droplet1 = sim.addDroplet(fluid1, dropletVolume, d_m1 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet2 = sim.addDroplet(fluid1, dropletVolume, d_m2 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet3 = sim.addDroplet(fluid1, dropletVolume, d_m3 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet4 = sim.addDroplet(fluid1, dropletVolume, d_m4 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet5 = sim.addDroplet(fluid1, dropletVolume, d_m5 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet6 = sim.addDroplet(fluid1, dropletVolume, d_m6 * flowRate / cWidth * cHeight, c1, 0.5);

    // check if chip is valid
    sim.checkChipValidity();

    // simulate
    droplet::SimulationResult result = sim.simulate();
    //std::cout << result.toJson(4) << std::endl;

    std::vector<int> expectedDropletHeaderPath = {1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36};
    std::vector<int> expectedDropletPayloadPath = {1, 3, 5, 6, 9, 11, 12, 15, 17, 18, 21, 23, 24, 27, 29, 30, 33, 35, 36};
    auto droplet0Path = result.getDropletPath(droplet0);
    auto droplet1Path = result.getDropletPath(droplet1);
    auto droplet2Path = result.getDropletPath(droplet2);
    auto droplet3Path = result.getDropletPath(droplet3);
    auto droplet4Path = result.getDropletPath(droplet4);
    auto droplet5Path = result.getDropletPath(droplet5);
    auto droplet6Path = result.getDropletPath(droplet6);

    auto check_path = [](droplet::DropletPath actualDropletPath, std::vector<int> expectedDropetPath) {
        int occuranceCount = 1;
        int idx = 0;
        for (auto dropletPath : actualDropletPath.positions) {
            for (auto channelId : dropletPath.channelIds) {
                ASSERT_NEAR(channelId, expectedDropetPath.at(idx), 1e-10);
                occuranceCount++;
                if (occuranceCount == 3) {
                    std::cout << channelId << " ";
                    occuranceCount = 0;
                    idx++;
                }
            }
        }
        //std::cout << std::endl;
    };
    check_path(droplet0Path, expectedDropletHeaderPath);
    check_path(droplet1Path, expectedDropletPayloadPath);
    check_path(droplet2Path, expectedDropletHeaderPath);
    check_path(droplet3Path, expectedDropletHeaderPath);
    check_path(droplet4Path, expectedDropletHeaderPath);
}

TEST(Droplet, RingNetworkE2) {
    //--- API ---
    droplet::Simulator sim;

    // flowRate pump
    auto flowRate = 3e-11;
    auto pump = sim.addFlowRatePump(-1, 0, flowRate);

    // droplet length
    auto lDroplet = 150e-6;

    // droplet distances
    double d_m1 = 881e-6;
    double d_m2 = 40678e-6;
    double d_m3 = 42185e-6;
    double d_m4 = 41475e-6;

    // channels
    auto cWidth = 100e-6;
    auto cHeight = 33e-6;
    auto bWidth = 200e-6;

    auto iLength = 2000e-6;

    // Node 1
    auto lIn1N1 = 469e-6;
    auto lIn2N1 = 531e-6;
    auto lByN1 = 1200e-6;
    auto lOut1N1 = 14531e-6;
    auto lOut2N1 = 24469e-6;
    // Node 2
    auto lIn1N2 = 743e-6;
    auto lIn2N2 = 757e-6;
    auto lByN2 = 1700e-6;
    auto lOut1N2 = 14257e-6;
    auto lOut2N2 = 24243e-6;
    // Node 3
    auto lIn1N3 = 1021e-6;
    auto lIn2N3 = 979e-6;
    auto lByN3 = 2200e-6;
    auto lOut1N3 = 13979e-6;
    auto lOut2N3 = 24021e-6;
    // Node 4
    auto lIn1N4 = 1300e-6;
    auto lIn2N4 = 1200e-6;
    auto lByN4 = 2700e-6;
    auto lOut1N4 = 13700e-6;
    auto lOut2N4 = 23800e-6;
    // Node 5
    auto lIn1N5 = 1580e-6;
    auto lIn2N5 = 1420e-6;
    auto lByN5 = 3200e-6;
    auto lOut1N5 = 13420e-6;
    auto lOut2N5 = 23580e-6;
    // Node 6
    auto lIn1N6 = 1862e-6;
    auto lIn2N6 = 1638e-6;
    auto lByN6 = 3700e-6;
    auto lOut1N6 = 13138e-6;
    auto lOut2N6 = 23362e-6;

    auto cLength = 500e-6;

    auto dLength1 = 881e-6;
    auto dLength2 = 40678e-6;
    auto dLength3 = 42185e-6;
    auto dLength4 = 41475e-6;

    auto c1 = sim.addChannel(0, 1, cHeight, cWidth, iLength);
    //Node 1
    auto c2_n1 = sim.addChannel(1, 2, cHeight, cWidth, lIn1N1);
    auto c3_n1 = sim.addChannel(1, 3, cHeight, cWidth, lIn2N1);
    auto c4_n1 = sim.addChannel(2, 4, cHeight, cWidth, lOut1N1);
    auto c5_n1 = sim.addChannel(3, 4, cHeight, cWidth, lOut2N1);
    auto c6_n1 = sim.addChannel(4, 21, cHeight, cWidth, cLength);
    auto c7_n1 = sim.addBypassChannel(2, 3, cHeight, bWidth, lByN1);
    //Node 2
    auto c2_n2 = sim.addChannel(21, 22, cHeight, cWidth, lIn1N2);
    auto c3_n2 = sim.addChannel(21, 23, cHeight, cWidth, lIn2N2);
    auto c4_n2 = sim.addChannel(22, 24, cHeight, cWidth, lOut1N2);
    auto c5_n2 = sim.addChannel(23, 24, cHeight, cWidth, lOut2N2);
    auto c6_n2 = sim.addChannel(24, 31, cHeight, cWidth, cLength);
    auto c7_n2 = sim.addBypassChannel(22, 23, cHeight, bWidth, lByN2);
    //Node 3
    auto c2_n3 = sim.addChannel(31, 32, cHeight, cWidth, lIn1N3);
    auto c3_n3 = sim.addChannel(31, 33, cHeight, cWidth, lIn2N3);
    auto c4_n3 = sim.addChannel(32, 34, cHeight, cWidth, lOut1N3);
    auto c5_n3 = sim.addChannel(33, 34, cHeight, cWidth, lOut2N3);
    auto c6_n3 = sim.addChannel(34, 41, cHeight, cWidth, cLength);
    auto c7_n3 = sim.addBypassChannel(32, 33, cHeight, bWidth, lByN3);
    //Node 4
    auto c2_n4 = sim.addChannel(41, 42, cHeight, cWidth, lIn1N4);
    auto c3_n4 = sim.addChannel(41, 43, cHeight, cWidth, lIn2N4);
    auto c4_n4 = sim.addChannel(42, 44, cHeight, cWidth, lOut1N4);
    auto c5_n4 = sim.addChannel(43, 44, cHeight, cWidth, lOut2N4);
    auto c6_n4 = sim.addChannel(44, 51, cHeight, cWidth, cLength);
    auto c7_n4 = sim.addBypassChannel(42, 43, cHeight, bWidth, lByN4);
    //Node 5
    auto c2_n5 = sim.addChannel(51, 52, cHeight, cWidth, lIn1N5);
    auto c3_n5 = sim.addChannel(51, 53, cHeight, cWidth, lIn2N5);
    auto c4_n5 = sim.addChannel(52, 54, cHeight, cWidth, lOut1N5);
    auto c5_n5 = sim.addChannel(53, 54, cHeight, cWidth, lOut2N5);
    auto c6_n5 = sim.addChannel(54, 61, cHeight, cWidth, cLength);
    auto c7_n5 = sim.addBypassChannel(52, 53, cHeight, bWidth, lByN5);
    //Node 6
    auto c2_n6 = sim.addChannel(61, 62, cHeight, cWidth, lIn1N6);
    auto c3_n6 = sim.addChannel(61, 63, cHeight, cWidth, lIn2N6);
    auto c4_n6 = sim.addChannel(62, 64, cHeight, cWidth, lOut1N6);
    auto c5_n6 = sim.addChannel(63, 64, cHeight, cWidth, lOut2N6);
    auto c6_n6 = sim.addChannel(64, -1, cHeight, cWidth, cLength);
    auto c7_n6 = sim.addBypassChannel(62, 63, cHeight, bWidth, lByN6);

    //--- sink ---
    sim.addSink(-1);
    //--- ground ---
    sim.addGround(-1);

    // fluids
    auto fluid0 = sim.addFluid(1e-3, 1e3, 0.0, 9e-10);
    auto fluid1 = sim.addFluid(4e-3, 1e3, 1.0, 9e-10);
    //--- continuousPhase ---
    sim.setContinuousPhase(fluid0);

    // droplets
    auto dropletVolume = lDroplet * cWidth * cHeight;
    auto droplet0 = sim.addDroplet(fluid1, dropletVolume, 0.0, c1, 0.5);
    auto droplet1 = sim.addDroplet(fluid1, dropletVolume, d_m1 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet2 = sim.addDroplet(fluid1, dropletVolume, d_m2 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet3 = sim.addDroplet(fluid1, dropletVolume, d_m3 * flowRate / cWidth * cHeight, c1, 0.5);
    auto droplet4 = sim.addDroplet(fluid1, dropletVolume, d_m4 * flowRate / cWidth * cHeight, c1, 0.5);

    // check if chip is valid
    sim.checkChipValidity();

    // simulate
    droplet::SimulationResult result = sim.simulate();
    //std::cout << result.toJson(4) << std::endl;

    std::vector<int> expectedDropletHeaderPath = {1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36};
    std::vector<int> expectedDropletPayloadPath = {1, 3, 5, 6, 8, 10, 12, 15, 17, 18, 21, 23, 24, 26, 28, 30, 33, 35, 36};
    auto droplet0Path = result.getDropletPath(droplet0);
    auto droplet1Path = result.getDropletPath(droplet1);
    auto droplet2Path = result.getDropletPath(droplet2);
    auto droplet3Path = result.getDropletPath(droplet3);
    auto droplet4Path = result.getDropletPath(droplet4);

    auto check_path = [](droplet::DropletPath actualDropletPath, std::vector<int> expectedDropetPath) {
        int occuranceCount = 1;
        int idx = 0;
        for (auto dropletPath : actualDropletPath.positions) {
            for (auto channelId : dropletPath.channelIds) {
                ASSERT_NEAR(channelId, expectedDropetPath.at(idx), 1e-10);
                occuranceCount++;
                if (occuranceCount == 3) {
                    std::cout << channelId << " ";
                    occuranceCount = 0;
                    idx++;
                }
            }
        }
        //std::cout << std::endl;
    };
    check_path(droplet0Path, expectedDropletHeaderPath);
    check_path(droplet1Path, expectedDropletPayloadPath);
    check_path(droplet2Path, expectedDropletHeaderPath);
    check_path(droplet3Path, expectedDropletHeaderPath);
    check_path(droplet4Path, expectedDropletHeaderPath);
}
