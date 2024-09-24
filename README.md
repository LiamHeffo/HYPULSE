The HYPULSE driver script along with an extensive library of examples using PITOT to simulate previous HYPULSE runs.

-- HYPULSE_driver.py contains the functions used to calculate the conditions driving the primary shock. 

-- There are three types of driver functions:
1. HYPULSE_driver(): Standard fixed volume driver operation. Takes in fill conditions and composition and returns the plateau state.
2. HYPULSE_driver_equivalentconditions(): Returns the equivalent sound speed and pressure for a passive burst fixed volume driver as per \cite{}.
3. HYPULSE_driver_exp_tuned(): Returns the effective sound speed and pressure based on the regressions done in Heffernan, L., James, C.M. and Jewell, J.S. (2024). HYPULSE Modelling Using Quasi-0D Analytical Framework. AIAA AVIATION FORUM AND ASCEND 2024. doi:https://doi.org/10.2514/6.2024-4566.

Each driver functions accepts the same inputs: the fill pressure, temperature, and composition for states 400 and 100.

The library of examples is split into two main directories:
1. The TsaiBakos_RST folder, which contains the following sub-directories:
   1. Experimental_Models: A list of HYPULSE RST conditins which use PITOT3 in its RST mode. That is, a primary shock speed is specified and the downstream conditions tuned to the primary shock speed.
   2. RST_equivalentconditions: A list of HYPULSE RST conditions that use the equivalent condition driver function to simulate the shock tunnel.
   3. Theoretical_Models: A list of HYPULSE RST conditions using the standard theoretical driver function.
2. TsaiBakos_SET, which contains:
    1. SET_CJstate: A list of HYPULSE SET conditions using the CJ state post-detonation as the driver gas.
   2. SET_equivalentconditions: A list of HYPULSE SET conditions using the equivalent condition driver function.
   3. SET_Plateau_conditions: A list of HYPULSE SET conditions using the standard plateau conditions post detonation as the driver gas.

The repository also containts the Parametric_Studies directory, which contains an example of a parametric study on the M8 RST condition, along with example plots. For more information on how to run parametric studies using PITOT3, please see the PITOT3 documentation on the GDTk website.

To run any of the files, you will need to add the HYPULSE_driver.py file to your PYTHOPATH. You will also need to have the GDTk installation {which version?} correctly installed.