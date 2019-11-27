**Analysis settings**

### Which level of cuts to apply:
cuts_level: 1

### Analysis category ( 2-tracks | 3-layers | 4-layers | 5-6-layers | all ):
analysis_category: all

### Secondary category ( Zmumu | Wmunu ):
secondary_category: 

### Should events after processing be saved on disk:
save_events:  1

### Do we need additional information stored in the friend trees
load_friend_tree: 0
load_hits:             0

verbosity_level: 0

### turn on/off different backgrounds, signals and data samples
do_QCD:         0
do_Zmm:         0
do_tops:          0
do_dibosons:   0
do_Wmv:         0
do_Zvv:           0

do_300_3:       0
do_300_10:     0
do_300_30:     1
do_500_10:     0
do_500_20:     0
do_650_10:     0
do_650_20:     0
do_800_10:     0
do_800_20:     0
do_1000_10:   0
do_1000_20:   0
do_tagger_signal_noPU:                0
do_tagger_background_noPU:      0
do_tagger_signal_withPU:             0
do_tagger_background_withPU:   0

do_chargino_300_1:   0
do_chargino_400_1:   0 
do_chargino_300_10: 0
do_chargino_500_10: 0
do_chargino_500_1:   0
do_chargino_700_10: 0
do_chargino_700_30: 0
do_chargino_800_10: 0

do_SR: 0
do_CR: 0

### For samples slit into several chunks, one can load all chunks or just a single one.
### Setting load_single_subpath: 1 at L1 means that for each chunk input files will be loaded
### from its corresponding directory, rather than from the first one.

load_single_subpath: 0
subpath_index: 0

### Select for which years to run:
load_2017: 1
load_2018: 1

**Printing & plotting options**

### Print total yields for background, data and signal before and after processing:
print_yields: 1

### Print detailed description of events passing selections:
print_background_details: 0
print_data_details: 0
print_signal_details: 0

### Draw standard (per event/track/jet/helix) plots:
draw_standard_plots:  0

### Draw per-layer plots (mainly dE/dx):
draw_per_layer_plots: 0

### Show legends in plots:
show_legends: 1

**Other analysis options**

### Limit number of events loaded (-1 means load all available)
max_N_events_background:  -1
max_N_events_signal: -1
max_N_events_data:  -1



### Luminosity (in fb^-1):

total_luminosity_2015: 3.81
total_luminosity_2016: 37.76
total_luminosity_2017: 41.37
total_luminosity_2018: 63.97

total_luminosity_run2: 146.91
