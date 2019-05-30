# PlanktonBS
Differential Equation-Based model for studying the seasonal succession of phytoplankton in the Eastern Bering Sea.

# Description
The model includes (1) four phytoplankton functional types (diatoms, flagellates, autotrophic dinoflagellates, and the coccolithoporid *Emiliania huxleyi*), (2) two zooplankton types (micrzooplankton and mesozooplankton), (3) three types of dissolved inorganic nutrients (nitrate, ammonium, and silicate), (4) detritus, (5) attached and free coccoliths, (6) dissolved inorganic carbon (DIC), and (7) total alkalinity, for a total of 14 differential equation. The model resolves fully the carbonate system by calculating (based on DIC and alkalinity): bicarbonate ion ocncentration, carbonate ion concentration, omega calcite, omega aragonite, and pH.

# How to run
The model has to be compiled with C++ from the [Gnu Compiler Collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection) using the command `g++` as follows:

```
g++ succession4new.cc routines.cc nrutil.cc -o a.out -Wno-deprecated
```

This creates the executable called `a.out`, which is run by typing `./a.out`.


