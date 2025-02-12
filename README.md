# Modular 1D-Simulator

[![Ubuntu CI](https://img.shields.io/github/actions/workflow/status/cda-tum/mmft-modular-1D-simulator/ubuntu.yml?label=Ubuntu&logo=ubuntu&style=flat-square)](https://github.com/cda-tum/mmft-modular-1D-simulator/actions/workflows/ubuntu.yml)
[![macOS CI](https://img.shields.io/github/actions/workflow/status/cda-tum/mmft-modular-1D-simulator/macos.yml?label=macOS&logo=apple&style=flat-square)](https://github.com/cda-tum/mmft-modular-1D-simulator/actions/workflows/macos.yml)
[![Nature](https://img.shields.io/static/v1?label=Nature&message=Modular-1D-Simulator&color=informational&style=flat-square)](https://doi.org/10.1038/s41598-024-77741-8)

<p align="center">
  <picture>
    <img src="https://www.cda.cit.tum.de/research/microfluidics/logo-microfluidics-toolkit.png" width="60%">
  </picture>
</p>

A 1D-Simulator for Microfluidic Biochips developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/) as part of the [Munich MicroFluidic Toolkit (MMFT)](https://www.cda.cit.tum.de/research/microfluidics/munich-microfluidics-toolkit/).
This 1D-Simulator exploits the 1D analysis model, which is especially suited for simulating designs before even the first prototype is fabricated and for design space explorations. Furthermore, it allows to simulate continious flow, instantaneous mixing, droplets and their respective paths, and membranes inside a Lab-on-a-Chip (LoC) with closed micro-channels. 
For more information about our work on Microfluidics, please visit https://www.cda.cit.tum.de/research/microfluidics/. 
More details about the implementation can be found in:

M. Emmerich, F. Costamoling, and R. Wille. [Modular and extendable 1D-simulation for microfluidic devices](https://doi.org/10.1038/s41598-024-77741-8). Scientific Reports, 2024.

G. Fink, F. Costamoling, and R. Wille. [MMFT Droplet Simulator: Efficient Simulation of Droplet-based Microfluidic Devices](https://doi.org/10.1016/j.simpa.2022.100440). Software Impacts, 2022.

If you have any questions, feel free to contact us via microfluidics.cda@xcit.tum.de or by creating an issue on GitHub. 

## System Requirements
The implementation should be compatible with any current C++ compiler supporting C++17 and a minimum CMake version 3.21. 

If you have doxygen installed, a detailed documentation can be build with the following command inside the build folder of the project: 
```bash 
make dropletDocumentation
```

The documentation will be written to the doc folder within your build folder.

## Usage
To use this library include the following code in your cmake file: 
```cmake
include(FetchContent)
FetchContent_Declare(
    droplet
    GIT_REPOSITORY https://github.com/cda-tum/mmft-droplet-simulator.git
    GIT_TAG master
)
FetchContent_MakeAvailable(droplet)

target_link_libraries(${TARGET} PRIVATE droplet)
```
and include the library API header in your project file:
```cpp
#include "droplet-simulator/Simulator.h"
```

## Step-by-Step Guide
This guide will walk you through the process of running the simulator, even if you have no prior experience with C++ programming.
In the future, we plan to introduce a GUI to simplify the process.

### Initial Setup
Before using the tool the first time make sure to complete the following steps:
1. **Install an IDE:** If you don't have one installed already, download and install an IDE that supports C++ development, such as CLion or Visual Studio Code.
2. **Install a C++ Compiler and CMake:** To compile C++ code, you'll need a compiler and build system. Here is how to set them up on different operating systems:
- Ubuntu(Linux): 
```bash
sudo apt update
sudo apt install build-essential cmake
```
- Windows:
```bash
choco install cmake
```
- MacOS:
```bash
xcode-select --install
brew install cmake
```
3. **Clone the Repository:** Use Git to clone the project repository to your local machine. You can do this using the Git command line or the GitHub interface (e.g., GitHub.com or GitHub Desktop).

### Step 1 - Build the executable
1. **Create a Build Folder and Compile:** Open your terminal and follow these command to create a build directory and compile the code
```c++
mkdir build
cd build
cmake ..
make
```

### Step 2 - Choose a Test Example and Run the Simulation
1. **Choose a Test Example:** Navigate to the `test` folder and select an example to use as a starting point for your simulation. We recommend runnning the code without any modifications initially to ensure everything is set up correctly. 
2. **Run the Test Simulation:** You can run any of the test cases using Google Test (also after you have adapted them). To run an example, such as the membrane test found in `tests/organ/Organ.tests.cpp`, use the following command in your terminal:
```c++
./dropletTest --gtest_filter=Organ.TwoOrgan
```
You can adapt the executable name and test case to match the simulation you want to run.

### Step 3 - Customize the Test File
1. **Select a Suitable Example:** Choose a test case that closely matches your desired simulation.
2. **Edit the Code:** Review the comments in the example code, which explain each line of code. Adjust the values or channel definitions as needed for your simulation.

If you create a new test file, don't forget to:
- Update the `CMakeLists.txt` file in your new test folder.
- Ensure that the new test is referenced in the `CMakeLists.txt`file located in the `tests`directory.

## Example
This small example shows how to create and simulate a small network for droplet simulations. 
```c++
#include <iostream>

#include "droplet-simulator/Results.h"
#include "droplet-simulator/Simulator.h"

int main(int argc, char const* argv[]) {
    std::cout << "--- Main ---" << std::endl;
    // create the simulator
    droplet::Simulator sim;

    std::cout << "--- flowRatePump ---" << std::endl;
    // flowRate pump
    auto flowRate = 3e-11;
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
    auto c4 = sim.addChannel(2, 4, cHeight, cWidth, cLength);
    // create channel from node 3 to node 4 with the given height, width, length
    auto c5 = sim.addChannel(3, 4, cHeight, cWidth, cLength);
    // create channel from node 4 to node -1 with the given height, width, length
    auto c6 = sim.addChannel(4, -1, cHeight, cWidth, cLength);
    std::cout << "--- sink ---" << std::endl;
    // define that node -1 is a sink
    sim.addSink(-1);
    std::cout << "--- ground ---" << std::endl;
    // define that node -1 is the ground node
    sim.addGround(-1);

    // fluids
    std::cout << "--- fluids ---" << std::endl;
    // add fluid with 1e-3 viscosity and 1e3 density
    auto fluid0 = sim.addFluid(1e-3, 1e3, 9e-10);
    // add fluid with 3e-3 viscosity and 1e3 density
    auto fluid1 = sim.addFluid(3e-3, 1e3, 9e-10);
    std::cout << "--- continuousPhase ---" << std::endl;
    // define which fluid is the continuous phase
    sim.setContinuousPhase(fluid0);

    // droplet
    std::cout << "--- droplet ---" << std::endl;
    auto dropletVolume = 1.5 * cWidth * cWidth * cHeight;
    // create a droplet of fluid1, with a given droplet volume, injected at injectionTime 0.0 in channel with the channelId c1 at relative injection position 0.5
    auto droplet0 = sim.addDroplet(fluid1, dropletVolume, 0.0, c1, 0.5);

    std::cout << "--- validity check chip ---" << std::endl;
    // check if chip is valid
    sim.checkChipValidity();

    std::cout << "--- simulate ---" << std::endl;
    // simulate the microfluidic network
    auto result = sim.simulate();

    // print the result
    std::cout << "--- result ---" << std::endl;
    std::cout << result.toJson(4) << std::endl;

    return 0;
}
```

More examples are defined in the tests subfolder.

