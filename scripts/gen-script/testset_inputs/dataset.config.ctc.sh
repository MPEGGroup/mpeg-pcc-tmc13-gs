#!/bin/bash
# List of all test sequences and their corresponding render resolutions.
# This file is sourced by the main script.
# The function ADD_SEQ(name, width, height) is already defined in the main script.

# NF and 1F sequences
ADD_SEQ bartender_semitracked   1920 1080
ADD_SEQ cinema_semitracked      1920 1080
ADD_SEQ breakfast_semitracked   1920 1080
ADD_SEQ breakfast_untracked     1920 1080
ADD_SEQ breakdance_untracked    1920 1080
ADD_SEQ manwithfruit_tracked    3840 2160

# Optional NF sequences
ADD_SEQ bartender_tracked       1920 1080
ADD_SEQ cinema_tracked          1920 1080
ADD_SEQ breakfast_tracked       1920 1080

# 1F sequences (only included if FLAG_1F == 1)
if [ "${FLAG_1F:-0}" -eq 1 ]; then
    ADD_SEQ trio                1920 1080
    ADD_SEQ makeup              1920 1080
    #ADD_SEQ musical             1920 1080
    ADD_SEQ lego_ferrari        4594 5514
    ADD_SEQ lego_bugatti        3852 2868
    ADD_SEQ cricket_player      4520 2540
    ADD_SEQ plant               2954 3968
    ADD_SEQ solo_tango_female   4534 2542
    ADD_SEQ solo_tango_male     4528 2544
    ADD_SEQ tango_duo           4522 2538
    ADD_SEQ tennis_player       4518 2540
    ADD_SEQ library             3793 2131
    ADD_SEQ flowerdance         2456 2054
    ADD_SEQ gymnast             2456 2054
fi
