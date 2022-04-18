# Apophis
## REBOUND/REBOUNDx Code for Close Approach Structure Analysis of Apophis

Tam Q. Do April 16, 2022

This repository houses the code used for my Honours Thesis project investigating the possible orbital evolution of the asteroid 99942 Apophis. All code is written in `python`, primarily using the N-body integration package `REBOUND` and the extra physics add-on package `REBOUNDx`. 

### apophis_template.py
This script is the main script used to run long-term integrations of Apophis' possible orbital trajectories using the LOV sampling method (with the possibility of evenly sampling its 7-dimensional uncertainty ellipsoid). For each simulation run, it prints useful information that can be used for mapping the close approach structure on a specific b-plane (more details on the output log can be found in [Output Log](#output-log)). For high-fidelity simulations, it considers 
* The gravitational perturbations from the Sun, Moon, planets in our Solar System, and 25 asteroids with the largest gravitational effect on Apophis
* The gravitational harmonics J2 and J4 due to the Earth's oblateness
* The radiation pressure force and Yarkovksy effect
* Post-Newtonian relativistic corrections for the Sun (can be changed to include the corrections for all bodies in the simulation, however, the run time of the script increases significantly and in most cases, the correction is only significant for the most massive body in the simulation)

The script is customizable and commented throughout for what variables can be changed for different simulation runs. The main customizations include:
* Length of the simulations in years (change variables `tsimend`, `dtout`)
* Which b-plane the close approach structure is examined on (change variables `tbfreeze_start`, `tbfreeze_end`)
* Velocity range sampled for LOV sample and sampling resolution (change variable `vpert`, `npert`)

Note that this script forces all simulations to start at the date January 21, 2022 since that was the date the solution and orbital parameters were calculated for by Davide Farnocchia on June 29, 2021. As such, we have made a one time pull of the orbital parameters for all bodies used in the simulations from the JPL ephemerides for the start date by uncommenting the code block in lines 82-117 (pulls info) and 136-139 (formats and prints info). We copy and paste the printed information into the code block in lines 152-187 to ensure all bodies start at the same start date. All bodies are integrated forwards in time after that point.

To run the file and save the output to a log file:
>: python3 apophis_template.py > logfile.txt

### comparison.py
This script is the secondary script used to compare the difference between nominal orbital trajectories under different considerations. That is, instead of sampling a range of possible trajectories, we simply run one simulation for one virtual Apophis asteroid. For each run, parts of the code can be commented and uncommented (this is well documented in the script itself) to run the simulation under the following scenarios:
* In addition to the planets, Sun and moon, only the gravitational perturbations from the 4-10 most massive asteroids are considered
* The post-Newtonian GR correction to the Sun is not considered
* The radiation pressure forces besides the Yarkovsky effect are not considered
* The gravitational harmonics due to Earth's oblateness (J2/J4) are not considered
* The Yarkovsky effect is modified at a time `tbfreeze`

To run the file and save the output to a log file:
>: python3 comparison.py > logfile.txt

### Output Log
The masses of the asteroids used (in solar masses) are printed once before starting to run simulations. For each simulation run, the following lines are printed:
* The perturbations pulled from the uncertainty ellipsoid (should be a list of 7 zeros if `SIG` is set to 0) and the perturbation-adjusted parameters for Apophis
* The LOV velocity perturbation and adjusted velocities
* The minimum distance between Apophis and Earth in the `tbfreeze_start` - `tbfreeze_end` time range and the time at which that distance is achieved (this is where the b-plane is frozen in)
* The distance and the zeta and xi coordinate on the frozen in b-plane for the location that Apophis passes through for that close approach
* The minimum distance after that close approach and the time at which that distance is achieved
* The incoming velocity of Apophis relative to Earth which is used to determine the b-plane orientation and evaluate the Earth's impact radius including gravitational focusing (`Vinf`)
* Correct Prediction which is just a sanity check of the Rebound collision predictor vs our distance/radius of impact calculations where 1 indicates agreement and 0 means we should probably check our calculations)
