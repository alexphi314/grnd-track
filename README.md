## Overview

This project generates a ground track for a given orbiting object and reports upcoming visible passes from a reference
(observation) location. It will also plot the ground track while the object is visible at the reference location.

## Installation

Python modules are defined in requirements.txt. Basemap is an extension of matplotlib that enables ground track plotting, 
requiring the GEOS library. For installation instructions, see https://matplotlib.org/basemap/users/installing.html.

## Running

In order to generate a ground track, run _initialize_track.py_. To see all options:
`python initialize_track -h`.

Initial conditions are required. To input Keplerian orbit elements, write a file in the following format (one per line):

epoch (UTC)
semi-major axis (m)
eccentricity
inclination (deg)
true anomaly (deg)
RAAN (deg)
argument of perigee (deg)

To input a TLE, supply a file with the TLE, or supply a NORAD satellite identifier, and the script will fetch the most
recent TLE from space-track.org. In order to fetch from space-track.org, you must supply user credentials as environment
variables SATCAT_USER and SATCAT_PASSWORD.

Example initial condition files are ex_ic.kep and ex_ic.tle.

## Getting observable passes

To view observable passes, you must supply the observation coordinates in a file. Use the following format:

<LAT (N)> <LON (E)>, in DMS
`39-45-43, -104-52-52` represents 39˚45'43'' N, -104˚52'52''E

Given that the initial conditions are sourced from TLEs, observable passes will only be generated for the next three days,
by default.

Example observation coordinates file is ex_rc.txt

## Viewing output

The reference ground track will be saved in the base project directory at _ground_track.png_. If specified, the visible 
pass plots will be found in the Plots/ folder.