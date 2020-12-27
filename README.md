# MESH-SED
Semi-distributed and physically-based sediment transport model for MESH. 

## Compilation
1. go into `src` 
2. type: `make`. It produces an executable `sa_mesh_sed` in `bin`
3. (optional) type `make clean` to delete `*.o` files and the executable `sa_mesh_sed`

## Execution
type: `./bin/sa_mesh_sed path/to/the/project` (e.g. `./bin/sa_mesh_sed test`). The`path/to/the/project` is a directoty that contain the following input files:

1. `MESH_sed_gridCell.ini` 
2. `MESH_sed_soilVegChan.ini` 
3. `MESH_sed_param.ini` 
4. `MESH_sed_reservoir.ini`

and the following output directories:

1. `output_field`
2. `output_ts` 

