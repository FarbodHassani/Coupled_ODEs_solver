#!/bin/bash
./one_orbit_5parameters "$1" "$2" "$3" -0.5 1.e-5 | perl ./plot5.pl | xmgrace -