/**
 * @file MembraneResistanceModels.h
 */

#pragma once

#include "IMembraneResistanceModel.h"
#include "Mixture.h"

namespace sim {

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 0
 * Based on Source: M. Ishahak, J. Hill, Q. Amin, L. Wubker, A. Hernandez, A. Mitrofanova, A. Sloan, A. Fornoni and A. Agarwal. "Modular Microphysiological System for Modeling of Biologic Barrier Function". In: Front. Bioeng. Biotechnol. 8:581163. (2020). doi: 10.3389/fbioe.2020.581163
 */
class MembraneResistanceModel0 : public sim::IMembraneResistanceModel {
  private:
    double getPoreResistance(arch::Membrane const* const membrane, Fluid const* const fluid) const;

  public:
    MembraneResistanceModel0();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
    double getPermeabilityParameter(const arch::Membrane* const membrane) const;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 1
 * Based on Source: J. J. VanDersarl, A. M. Xu, and N. A. Melosh. “Rapid spatial and temporal controlled signal delivery over large cell culture areas”. In: Lab Chip 11 (18 2011), pp. 3057–3063. doi: 10.1039/C1LC20311H. url: http://dx.doi.org/10. 1039/C1LC20311H.
 */
class MembraneResistanceModel1 : public sim::IMembraneResistanceModel {
  private:
    double diffusionCoefficient;

    double getPoreExitResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;
    double getPoreResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;
    double getPoreDensityDependentResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel1();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 2
 * Based on Source: J.L. Snyder, A. Clark, D.Z. Fang, T.R. Gaborski, C.C. Striemer, P.M. Fauchet, J.L. McGrath, "An experimental and theoretical analysis of molecular separations by diffusion through ultrathin nanoporous membranes". Journal of Membrane Science Volume 369, Issues 1–2, 2011, Pages 119-129. ISSN 0376-7388,https://doi.org/10.1016/j.memsci.2010.11.056.
 */
class MembraneResistanceModel2 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double area, double freeDiffusionCoefficient) const;

    double getTransmembraneResistance(arch::Membrane const* const membrane, double diffusionCoefficient, double area) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel2();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 3
 * Assumes N disk-like absorbers on the surface of a sphere
 * Based on Source: H.C. Berg. Random Walks in Biology. Princeton paperbacks. Princeton University Press, 1993. isbn: 9780691000640. url: https://books.google.at/books?id=DjdgXGLoJY8C
 */
class MembraneResistanceModel3 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient, double area) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel3();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 4
 * Based on Source: C.-W. Ho, A. Ruehli, and P. Brennan. “The modified nodal approach to network analysis”. In: IEEE Transactions on Circuits and Systems 22.6 (1975), pp. 504–509. doi: 10.1109/TCS.1975.1084079.
 */
class MembraneResistanceModel4 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel4();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 5
 * Adapted from Membrane Resistance Model 4 to fit the calculations in Appendix 1
 * Based on Source: C.-W. Ho, A. Ruehli, and P. Brennan. “The modified nodal approach to network analysis”. In: IEEE Transactions on Circuits and Systems 22.6 (1975), pp. 504–509. doi: 10.1109/TCS.1975.1084079.
 */
class MembraneResistanceModel5 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel5();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 6
 * Based on Source: K. Ronaldson-Bouchard et al. “A multi-organ chip with matured tissue niches linked by vascular flow”. In: Nature Biomedical Engineering 6.4 (2022), pp. 351–371. doi: 10.1038/s41551-022-00882-6. url: https://doi.org/10.1038/s41551-022-00882-6 (cit. on p. 35).
 */
class MembraneResistanceModel6 : public sim::IMembraneResistanceModel {
  private:
    double diffusionFactorMembrane = 1.58e-11;  // MultiOrgan

    //double diffusionFactorMembrane = 4.4e-11;  // TwoOrgan

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel6();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 7
 * Based on Source: H.C. Berg. Random Walks in Biology. Princeton paperbacks. Princeton University Press, 1993. isbn: 9780691000640. url: https://books.google.at/books?id=DjdgXGLoJY8C
 */
class MembraneResistanceModel7 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel7();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 8
 * Assumes N circular pores in a planar barrier, with the distance as membrane height
 * Based on Source: H.C. Berg. Random Walks in Biology. Princeton paperbacks. Princeton University Press, 1993. isbn: 9780691000640. url: https://books.google.at/books?id=DjdgXGLoJY8C
 */
class MembraneResistanceModel8 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel8();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

/**
 * @brief Class that defines the functionality of the 1D membrane resistance model.
 * Membrane Resistance Model 9
 * Assumes N circular pores in a planar barrier, with the distance as membrane height, and the diffusion resistance equal to a disk-like absorber
 * Based on Source: H.C. Berg. Random Walks in Biology. Princeton paperbacks. Princeton University Press, 1993. isbn: 9780691000640. url: https://books.google.at/books?id=DjdgXGLoJY8C
 */
class MembraneResistanceModel9 : public sim::IMembraneResistanceModel {
  private:
    double getPoreDiscoveryResistance(arch::Membrane const* const membrane, double diffusionCoefficient) const;

    double getPorePassageResistance(arch::Membrane const* const membrane, double area, double diffusionCoefficient) const;

  public:
    /**
     * @brief Instantiate the resistance model.
     */
    MembraneResistanceModel9();

    double getMembraneResistance(arch::Membrane const* const membrane, Fluid const* const fluid, double area) const override;
};

}  // namespace sim