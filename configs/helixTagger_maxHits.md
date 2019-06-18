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

**Pion's parameters**

### Minimum and maximum allowed momentum:
min_px: 50
min_py: 50
min_pz: 10

max_px: 250
max_py: 250
max_pz: 1000

**Chargino's track parameters**

### Maximum track's pseudorapidity:
max_eta:  1.2

### Number of tracker hits left by the chargino (only for generating random track, otherwise it's read from the event):
n_track_hits: 2

**Fitter parameters**

### Determines how far points can be from helix to be assigned to it (in mm):
helix_thickness:  0.0

### Maximum distance in 3D to consider two hits as overlapping ones (in mm):
double_hit_max_distance:            16

### Constraints on seeds parameters:
seed_max_chi2:                            3

seed_middle_hit_min_delta_phi:   -0.1
seed_middle_hit_max_delta_phi:  0.1
seed_middle_hit_max_delta_z:     60

seed_last_hit_min_delta_phi:       -0.8
seed_last_hit_max_delta_phi:       0.1
seed_last_hit_max_delta_z:          180

### Constrains on pion track parameters:
track_max_chi2:                             0.001

next_point_min_delta_phi:              -0.5
next_point_max_delta_phi:             0.5
next_point_max_delta_z:                220
next_point_max_delta_xy:              20

track_min_n_points:                       3

### Max number of different points and min number of points to merge two helices:
merging_max_different_point:         2
candidate_min_n_points:                3
merge_at_turn_back:                      0
merge_final_helices:                        1

### Max number of missing hits (total and in a row):
max_n_missing_hits:                       1
max_n_missing_hits_in_raw:           1

### Asymmetric hits constraints:
do_asymmetric_constraints:           1

### Turn on/off turning back helices
allow_turning_back:                         1

**Benchmark parameters**

### Conditinos to accept fitted helix as a proper solution (each parameter has to be within given tolerance).
### Position tolerance in mm, momentum tolerance in MeV:
tolerance_x:  10
tolerance_y:  10
tolerance_z:  10
tolerance_px:  10
tolerance_py:  10
tolerance_pz:  10
