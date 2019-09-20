**Analysis settings**

### Which level of cuts to apply:
cuts_level: 1

### Analysis category ( 2-tracks | 3-layers | 4-layers | all ):
analysis_category:  all

scan_MET_binning: 0
do_MET_binning: 0

### Should events after processing be saved on disk:
save_events:  1

### Do we need additional information stored in the friend trees
load_friend_tree: 0
load_hits:             0

verbosity_level: 1

### turn on/off different backgrounds, signals and data samples
do_QCD:         1
do_Zmm:         1
do_tops:          1
do_dibosons:   1
do_Wmv:         1
do_Zvv:           1

do_300_3:       0
do_300_10:     0
do_300_30:     0
do_500_10:     0
do_500_20:     0
do_650_10:     0
do_650_20:     0
do_800_10:     0
do_800_20:     0
do_1000_10:   0
do_1000_20:   0

do_2017:         0

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
show_legends: 0

**Other analysis options**

### Limit number of events loaded (-1 means load all available)
max_N_events_background:  -1
max_N_events_signal: -1
max_N_events_data:  -1



### Luminosity (in fb^-1) [2015: 3.81, 2016: 37.76, 2017: 41.37, 2018: 63.97, total Run 2: 146.91]:
total_luminosity: 146.91
