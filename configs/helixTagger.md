**General settings**

### Path of the output file:
output_path: helixFittingResults/tests.root

verbosity_level:  0

### Limit number of events loaded (-1 means load all available)
max_N_events_signal:  400

**Fitter parameters**

### Maximum distance in 3D to consider two hits as overlapping ones (in mm):
double_hit_max_distance:            20.0

### Constraints on seeds parameters:
seed_max_chi2:                             0.08

seed_middle_hit_min_delta_phi:   -0.6
seed_middle_hit_max_delta_phi:   0.9
seed_middle_hit_max_delta_z:     150

seed_last_hit_min_delta_phi:       -0.6
seed_last_hit_max_delta_phi:       0.7
seed_last_hit_max_delta_z:          250

### Constrains on pion track parameters:
track_max_chi2:                             0.01

next_point_min_delta_phi:             -0.6
next_point_max_delta_phi:             0.5
next_point_max_delta_z:                300
next_point_max_delta_xy:              250
next_point_max_delta_t:                0.6

track_min_n_points:                       3
track_min_n_layers:                       2

### Use distance to helix only when it passed through at least N layers:
min_layers_for_delta_xy:                5

### Max number of different points and min number of points to merge two helices:
merging_max_different_point:         2
candidate_min_n_points:                3
merge_at_turn_back:                      0
merge_final_helices:                        1

### Max number of missing hits (total and in a row):
max_n_missing_hits:                       1
max_n_missing_hits_in_raw:           1

### Asymmetric hits constraints:
do_asymmetric_constraints:           0

### Turn on/off turning back helices
allow_turning_back:                         1

### Reject seed if it causes starting values outside of limits:
require_good_starting_values:        1

### R(t) and s(t) functions:
exp_radius_function:        0
exp_slope_function:         0

### Let pion helix start one layer before/after the chargino track or have opposite charge
allow_one_less_layer: 0
allow_one_more_layer: 1
check_opposite_charge_below_Nlayers: 5

### Limits of helix parameters:
start_R0: 320
min_R0: 50
max_R0: 1000
min_Rslope:  0
max_Rslope: 10000

min_S0: -10000
max_S0: 10000
min_Sslope: -1000
max_Sslope: 0

min_X0: -5000
max_X0: 5000
min_Y0: -5000
max_Y0: 5000
min_Z0: -5000
max_Z0: 5000

**Helix tagger options**

### turn on/off different backgrounds, signals and data samples
do_QCD:         0
do_Zmm:         0
do_tops:          0
do_dibosons:   0
do_Wmv:         0
do_Zvv:           0

do_300_3:       0
do_300_10:     1
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

### Do we need additional information stored in the friend trees
load_friend_tree: 1


