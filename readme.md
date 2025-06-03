# Simulate time-dependent head in a semi-infinite leaky aquifer adjacent to surface water.

C.H. van Immerzeel

2/6/2025

Simulate time-dependent head in a semi-infinite leaky aquifer adjacent to 
open water (for example a river), where the boundary condition changes from 
t = 0 according to a specified course.

The equation used describes one-dimensional flow for a situation with 
spatially constant kD, c and S.

The response to a sequence of stress values (h) may be studied (simulation option).

![Formulas](https://github.com/user-attachments/assets/ae4afc5f-2d70-47a8-adbf-a62da07359bf)


## Link to the app
<https://sweco.shinyapps.io/Brug1DLeaky/>

## Source code
R-source code of the app:

<https://github.com/KeesVanImmerzeel/Brug1DLeaky>

## Simulation
A sequence of stress values (a) may be uploaded. For this purpose, prepare a spreadsheet with two columns ('Time' and 'h').  


## Numerical stability
Be aware that exotic input values may lead to numerical instability. In that case, no results or plots are presented. Instead, an error message appears.


## References
- Analytical solutions of geohydrological problems. G.A. Bruggeman p. 66 (Formula 123.32).


