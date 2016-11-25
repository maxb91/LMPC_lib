# barc_lib

This repo contains functions and types that are used by both the Barc simulation framework and the "real" ROS framework.
It contains basic BARC-specific types (related to the model and to MPC functions) and functions for their initialization. There are also functions to simulate the car's dynamics and to control the car along a racetrack (using Learning Model Predictive Control).
The library should be structured in folders according to its purposes.

History:
========
1. Introduced general functions and classes that can be used in simulation and on the BARC
2. Added specific functions for Learning Model Predictive Control (LMPC) and system identification
3. Added the option for steering delay (in the MPC formulation, approx. 2-step-delay on the BARC system (~0.2 seconds))
4. Added branch for system ID with a different sampling frequency than the MPC formulation (50 Hz sysID, 10 Hz MPC)