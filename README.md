Simple OWC MATLAB Model

This MATLAB script estimates the power output and efficiency of an Oscillating Water Column (OWC) wave-energy device using a simplified linear model.

Part 1 – Power and Efficiency Calculation

Given user-defined parameters (wave height H, period T, chamber and turbine areas Ac′, A′ₜ, and power coefficient Cₚ), the script:

Computes the vertical water and air velocities.

Calculates instantaneous and mean turbine power.

Compares mean turbine power (W/m) to theoretical wave power per meter to obtain efficiency.

Part 2 – Parametric Analysis

Tests how efficiency changes when varying:

Wave height H

Wave period T

Area ratio Ac′/At′

and plots η versus each parameter.
