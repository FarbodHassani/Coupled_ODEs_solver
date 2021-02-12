#!/bin/bash
./one_orbit_no_taustar "$1" "$2" "$3" | perl ./plot_one.pl | xmgrace -