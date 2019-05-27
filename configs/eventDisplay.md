**Event display options**

show_geometry_pixel:  0
show_geometry_strip:  0
show_geometry_ecal:  0
show_geometry_hcal:  0

draw_tracker_clusters:      0
draw_met:                         0
draw_jets:                          1
draw_pion_simhits:           0
draw_pion_clusters:          0
draw_chargino_simhits:    1

**Input options**

### turn on/off different backgrounds, signals and data samples
do_QCD:         0
do_Zmm:         0
do_tops:          0
do_dibosons:  0
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

### Determines how far points can be from helix to be assigned to it (in mm):
helix_thickness:  0.0

### Determines how far points can be from circle to be assigned to it (in mm):
circle_thickness: 20.0

### Constraints on seeds parameters:
seed_max_chi2:                            1e-6
seed_middle_hit_max_delta_phi:  1.0
seed_middle_hit_max_delta_z:     100
seed_last_hit_max_delta_phi:        1.0
seed_last_hit_max_delta_z:           100

### Constrains on pion track parameters:
track_max_chi2:                             1e-6
next_point_max_delta_phi:             1.0
next_point_max_delta_z:                200
track_min_n_points:                       5

### Max number of different points to merge two helices:
merging_max_different_point:        2

### Min number of points for a candidate to pass to the merging step:
candidate_min_n_points:                5

**Benchmark parameters**

### Conditinos to accept fitted helix as a proper solution (each parameter has to be within given tolerance).
### Position tolerance in mm, momentum tolerance in MeV:
tolerance_x:  10
tolerance_y:  10
tolerance_z:  10
tolerance_px:  10
tolerance_py:  10
tolerance_pz:  10

