**Input options**

### Which level of cuts to apply:
cuts_level: 1

### Analysis category ( 2-tracks | 3-layers | 4-layers | 5-6-layers | all ):
analysis_category: all

### Secondary category ( Zmumu | Wmunu | LowMET | none ):
secondary_category:  none

### turn on/off different backgrounds, signals and data samples
do_QCD:         0
do_Zmm:         0
do_tops:          0
do_dibosons:  0
do_Wmv:         0
do_Zvv:           0

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

do_tagger_signal_noPU:                0
do_tagger_background_noPU:      0
do_tagger_signal_withPU:             0
do_tagger_background_withPU:   0

do_chargino_300_1: 0
do_chargino_400_1: 0
do_chargino_500_1: 0
do_chargino_500_10: 1
do_chargino_700_10: 0
do_chargino_800_10: 0
do_chargino_700_30: 0

do_SR:         0
do_CR:         0

### Select for which years to run:
load_2017: 1
load_2018: 0

### Do we need additional information stored in the friend trees
load_friend_tree: 1
load_hits:             1

**General settings**

### Maximum distance in 3D to consider two hits as overlapping ones (in mm):
double_hit_max_distance: 20.0

### Limit number of events loaded (-1 means load all available)
max_N_events_signal:  -1

verbosity_level: 2
