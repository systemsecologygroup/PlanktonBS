# PlanktonBS
Differential Equation-Based model for studying the seasonal succession of phytoplankton in the Eastern Bering Sea.

# Location
The model is applied in a point location (M2) of the Southeastern Bering Sea, as indicated in the map below.

<p align="center">
  <img src="map.png" width="500">
</p>


# Model description
The model includes (1) four phytoplankton functional types (diatoms, flagellates, autotrophic dinoflagellates, and the coccolithoporid *Emiliania huxleyi*), (2) two zooplankton types (micrzooplankton and mesozooplankton), (3) three types of dissolved inorganic nutrients (nitrate, ammonium, and silicate), (4) detritus, (5) attached and free coccoliths, (6) dissolved inorganic carbon (DIC), and (7) total alkalinity, for a total of 14 differential equation. The model resolves the full carbonate system by calculating (based on DIC and alkalinity): bicarbonate ion concentration, carbonate ion concentration, omega calcite, omega aragonite, and pH. A simplified model schematic is shown below.

<p align="center">
  <img src="schematic.png" width="500">
</p>


# How to run the model
The model has to be compiled with C++ from the [Gnu Compiler Collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection) using the command `g++` as follows:

```
     g++ succession4new.cc routines.cc nrutil.cc -o a.out -Wno-deprecated
```

This creates the executable called `a.out`, which is run by typing `./a.out`. The option `-Wno-deprecated` avoid getting warnings about the usage of deprecated features.

The model requires input files (forcing environmental functions, including Mixed Layer Depth, Seas Surface Temperature, Wind Speed, and Salinity), which have to be stored in a subdirectory called `./input`.

Header files (`param.h` and `nrutil.h` have to be present in the current directory).

Crucial model parameters are:.

```c++
     # define NEQ 14         // number of ordinary differential equations
     # define Y 9            // number of years for which run the model (0 is one year cycle)
     # define IGNY 0         // number of years required by the model to reach equilibrium (spin-up)
     # define HOFY 4320      // hour of the year to consider for poincare' sections
```

Results are saved in a subdirectory called `results`.

# Related publication
This model was used in the following papers:

- A. Merico, T. Tyrrell, E. J. Lessard, T. Oguz, P. J. Stabeno, S. I. Zeeman, and T. E. Whitledge. [Modelling phytoplankton succession on the Bering Sea shelf: role of climate influences and trophic interactions in generating *Emiliania huxleyi* blooms 1997-2000](https://www.sciencedirect.com/science/article/pii/S0967063704001475). *Deep-Sea Research I*, **51**:1803–1826, 2004.

- A. Merico, T. Tyrrell, and T. Cokacar. [Is there any relationship between phytoplankton seasonal dynamics and the carbonate system?](https://www.sciencedirect.com/science/article/pii/S0924796305001892) *Journal of Marine Systems*, **59**:120–142, 2006.
