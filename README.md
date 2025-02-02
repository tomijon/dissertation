# Dissertation Project
## Investigating the performance of generating terrain on the CPU vs on the GPU.

This project investigates the performance difference of procedurally generating terrain in a real time system on the CPU vs on the GPU. It implements perlin noise and value noise,
and generates the terrain mesh using the marching cubes algorithm. The study concluded that the GPU is much faster at generating noise.
However, some implementation details meant the GPU version was slower at drawing. For example, there were a lot of empty triangles being drawn on the GPU version as I did not know about
geometry shaders. This basically tanks the FPS which is unfortunate but a possible improvement.
