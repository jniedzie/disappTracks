**General settings**
### Add hits coming from random pion to the collection of all hits (0 - false, 1 - true):
inject_pion_hits: 1

### Number of events to generate and process:
n_tests:  20

### Path of the output file:
output_path: helixFittingResults/tests.root

### Number of noise hits to generate in the pixel barrel:
n_noise_hits: 0

### Number of available tracker layers (this is used to calculate hits on silicon and generate noise):
n_tracker_layers: 10

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
helix_thickness:  0.3

### Determines how far points can be from circle to be assigned to it (in mm):
circle_thickness: 0.4

### Determines how far points can be from each other to be counted as the same z-line (in mm),
### when assigining circles (probably this needs to be greater than circle_thickness):
lines_tolerance_for_circles: 0.5

### Determines how far points can be from each other to be counted as the same z-line (in mm), 
### when calculating number of regular points:
lines_tolerance_for_regularity: 10.0

### Precision of pz scanning:
step_pz:  0.5

### Determines how far points can be from perfectly regular position in Z to be counted as a regular point (in mm):
z_regularity_tolerance: 1.0

### Minimum number of hits for each line along Z axis (approximately equivalent to minimum number of helix cycles).
### The higher this number, the faster it gets, but less cases will be successfully fitted:
min_n_points_along_z: 2

**Benchmark parameters**

### Conditinos to accept fitted helix as a proper solution (each parameter has to be within given tolerance).
### Position tolerance in mm, momentum tolerance in MeV:
tolerance_x:  10
tolerance_y:  10
tolerance_z:  10
tolerance_px:  10
tolerance_py:  10
tolerance_pz:  10