# Read Me

## Summary of Cases
This repository includes parametric sweep results conducted in the studies described in [published work] ...

Each folder describes one reflectivity distribution applied to a square channel honeycomb receiver:

1. Base Case: a uniform reflectivity receiver

2. Linear Distribution (Diffuse): a linear distribution case where the reflectivity increases from its lowest limit in the front of the receiver and starts ramping up at a slope $m$ after an entry length $L_e$

3. Linear Distribution (Specular): similar to the diffuse linear case, but the slope is large enough that the distribution is almost stepped instead of linear and the reflective region is specularly reflective

4. Wall Distribution: high reflectivtiy is applied to a pair of parallel walls, while the other pair of channel walls has the low reflectivtiy limit applied to them. 

## Repository Structure/Contents

The folder for each case (except the base case*) contains 3 sub-folders:

1. Stage 1: where optical properties (one reflectivity distribution variable $\bar(\varepsilon)$ and channel/pore radius $R_{ch}$) are varied at a fixed $L/D$ ratio of 25 and a $P/\dot{m}$ ratio of 700 kJ/kg. This folder will, also, contain the axial boundary heat source profiles for every combination of optical properties applicable to that design case.

2. Stage 2: where thermal properties ($L/D$ ratio) is varied at the ($\bar(\varepsilon)$, $R_{ch}$) pair optimized from stage 1 at $P/\dot{m}$ ratio of 700 kJ/kg 

3. Performance Chart: where the operating condition or $P/\dot{m}$ ratio is varied from 100 to 700 kJ/kg for the receiver with the optimized $\bar(\varepsilon)$, $R_{ch}$ and $L/D$ from stages 1 and 2. This folder will also contain lists of the extracted volumetric effct definitions and the notebook they were generated with.

Each of these folders, also, contains the list of parameters the model was setup with, the parametric sweep levels, temperature profile data, a `Global Powers.txt` file that contains the energy balance components extracted from COMSOL and the notebook/s where these text files are processed and visualized.

* The stage 1 folder for the base case includes two subdirectories. One includes the results of a study on the independent effect of porosity and pore radius at a fixed emissivity of 0.8, while the other folder contains the study of $\bar(\varepsilon)$ vs. pore radius like all the other receiver cases


## Usage & Dependencies

Apart from the text files, which could be downloaded and processed independently, one can build on the provided Pluto notebooks with only the following requirements:

- `julia >= 1.7`
- `Pluto.jl >= v"0.17.5"`
- `Plots.jl >= v"1.38.8"`

Simply download [julia](https://julialang.org/), launch a julia terminal, and run: 

```
import Pluto
Pluto.run()
```

and paste the path of the notebook you want to run in the provided location.
