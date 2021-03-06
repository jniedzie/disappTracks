**Event display options**

show_geometry_pixel:  0
show_geometry_strip:  0
show_geometry_ecal:  0
show_geometry_hcal:  0

draw_met:                         0
draw_jets:                          1

draw_tracker_clusters:      0
draw_pion_simhits:           0
draw_pion_clusters:          1
draw_chargino_simhits:    0

draw_true_helices:                    1
draw_fitted_helices:                  1
draw_fitted_helices_clusters:    1

fit_helices:                         1
fit_pion_clusters_only:       1
fit_noise_clusters_only:     0
include_endcaps:              1

**Input options**

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

do_chargino_500_10: 1

load_2017: 1

do_SR:         0
do_CR:         0

### Do we need additional information stored in the friend trees
load_friend_tree: 1
load_hits:             1

**General settings**
### Add hits coming from random pion to the collection of all hits (0 - false, 1 - true):
inject_pion_hits: 0

### Number of events to generate and process:
n_tests:  20

### Path of the output file:
output_path: helixFittingResults/tests.root

### Number of noise hits to generate in the pixel barrel:
n_noise_hits: 0

### Number of available tracker layers (this is used to calculate hits on silicon and generate noise):
n_tracker_layers: 4

### Limit number of events loaded (-1 means load all available)
max_N_events_signal:  -1

verbosity_level: 3

**Pion's parameters**

### Minimum and maximum allowed momentum:
min_px: 0
min_py: 0
min_pz: 100

max_px: 1000
max_py: 1000
max_pz: 1000

**Chargino's track parameters**

### Maximum track's pseudorapidity:
max_eta:  1.2

### Number of tracker hits left by the chargino (only for generating random track, otherwise it's read from the event):
n_track_hits: 2

**Fitter parameters**

### Maximum distance in 3D to consider two hits as overlapping ones (in mm):
double_hit_max_distance:            20.0

### Constraints on seeds parameters:
seed_max_chi2:                           0.01

seed_middle_hit_min_delta_phi:   -0.4
seed_middle_hit_max_delta_phi:  0.8
seed_middle_hit_max_delta_z:     150

seed_last_hit_min_delta_phi:       -0.4
seed_last_hit_max_delta_phi:       0.8
seed_last_hit_max_delta_z:          150

### Constrains on pion track parameters:
track_max_chi2:                             0.005

next_point_min_delta_phi:             -0.5
next_point_max_delta_phi:            1.5
next_point_max_delta_z:                500
next_point_max_delta_xy:              500
next_point_max_delta_t:                0.7

track_min_n_points:                       2
track_min_n_layers:                       3

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

### Let pion helix start one layer before/after the chargino track
allow_one_less_layer: 0
allow_one_more_layer: 1
check_opposite_charge_below_Nlayers: 1

### Limits of helix parameters:
start_R0: 320
min_R0: 50
max_R0: 1000
min_Rslope:  0
max_Rslope: 10000

min_S0: -10000
max_S0: 10000
min_Sslope: -150
max_Sslope: 0

min_X0: -5000
max_X0: 5000
min_Y0: -5000
max_Y0: 5000
min_Z0: -5000
max_Z0: 5000

**Benchmark parameters**

### Conditinos to accept fitted helix as a proper solution (each parameter has to be within given tolerance).
### Position tolerance in mm, momentum tolerance in MeV:
tolerance_x:  10
tolerance_y:  10
tolerance_z:  10
tolerance_px:  10
tolerance_py:  10
tolerance_pz:  10

