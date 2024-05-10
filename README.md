# Water Surface Wavelets

C++ implementation of the Water Surface Wavelets paper by Nvidia. 

Here are what our final results look like:

<img src = "walkthroughs/windsim.gif" width ="400" /> <img src = "walkthroughs/boat.gif" width ="400" />

These demonstrate our implementations of the core water surface simulation algorithm (advection + height field), profile buffer, diffusion, dissipation, procedurally generated terrain, solid-fluid coupling, and wind simulation that changes direction and intensity periodically.

## Controls

- Mouse/WASD for camera movement
- C to toggle orbit mode
- Arrow keys to control boat
- T to toggle wind

## Team Responsibilities

- All of us — Understanding paper math, defining amplitude interfaces, breaking down interpolation for advection step
- Max — Bilinear Interpolation, Diffusion Terms (Spatial and Angular), Terrain Generation, Amplitude initialization, framework for Boundary Interactions
- Dustin — Profile buffer, bilinear interpolation contributions, debugging advection step, and other parts of the code 
- Mohammed — Integrate system-solver for falling objects, solid-fluid coupling, boat movement inducing waves, adding wind sim w/ dissipation, a bit of reflective grid boundary conditions
- Kelvin — Basic interpolation, wave shading and rendering, skybox
