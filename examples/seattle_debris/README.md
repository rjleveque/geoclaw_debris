# seattle debris

the current code needs `SeattleBuildings/Seattle_Inventory.csv` and `SeattleBuildings/seattle_19arcsec_navd88_wgs84.tif`. the csv comes from `InventoryBRAILS_Seattle2.ipynb`. the tif must already be in EPSG:4326, in meters, and cover the domain in `setrun.py`. this version does not read the NOAA mhw `.nc` file.

with clawpack, its geoclaw source and makefiles, a fortran compiler, and the python dependencies installed, run from the repository root:

```sh
export CLAW=/path/to/clawpack
cd examples/seattle_debris
make flow
```

`make flow` calls `prepare_topo.py` to make the terrain, `make_qinit.py` to make the initial surge, geoclaw to run the model, and `combine_fgout_nc.py` to make the flow `.nc` file. `prepare_topo.py` writes `terrain_provenance.json` beside itself; `make flow` copies it into `_output` after the model run.

the defaults are a 12 m grid, 5 m surge, and 600 s run. the settings used by `seattle_tracer_flow.ipynb` are:

```sh
GEOCLAW_DX_METERS=6 GEOCLAW_SURGE_HEIGHT_M=12 GEOCLAW_TFINAL_SECONDS=900 make flow
```
