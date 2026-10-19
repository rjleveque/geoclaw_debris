# Seattle debris experiment

This experiment makes a GeoClaw flow field for the Seattle waterfront, then uses
it in the tracer and ship notebooks. The model terrain comes from [here](https://www.ngdc.noaa.gov/thredds/catalog/tiles/nthmp/tiled_19as/catalog.html?dataset=tilesDatasetScan/nthmp/tiled_19as/puget_sound_mhw_19_n47x75_w122x50_2023v1.nc).

## Inputs

1. Run `SeattleBuildings/InventoryBRAILS_Seattle2.ipynb` to create
   `SeattleBuildings/Seattle_Inventory.csv`. It must have a `geometry` column
   containing WKT building footprints in longitude and latitude.
2. `prepare_topo.py` uses
   [aforementioned data](https://www.ngdc.noaa.gov/thredds/catalog/tiles/nthmp/tiled_19as/catalog.html?dataset=tilesDatasetScan/nthmp/tiled_19as/puget_sound_mhw_19_n47x75_w122x50_2023v1.nc).
   Save `puget_sound_mhw_19_n47x75_w122x50_2023v1.nc` in
   `examples/seattle_debris`.

## Run

Install the Python packages in the repository `requirements.txt`, Clawpack
with its GeoClaw source and makefiles, and a Fortran compiler.

From the repo root:

```sh
export CLAW=/path/to/clawpack
cd examples/seattle_debris
make flow
```

`make flow` runs:

1. `prepare_topo.py`, which writes `seattle_19arcsec_mhw.tt3` and `seattle_buildings_mhw.tt3`
2. `make_qinit.py`, which creates the initial surge
3. GeoClaw, which writes frames to `_output`
4. `combine_fgout_nc.py`, which writes `pike_place_surge*_fgout.nc` for notebook use

The default run has a 12m grid, 5m surge, and 600s duration. The
settings used in `seattle_tracer_flow.ipynb` are

```sh
GEOCLAW_DX_METERS=6 GEOCLAW_SURGE_HEIGHT_M=12 GEOCLAW_TFINAL_SECONDS=900 make flow
```

Then you can run `seattle_tracer_flow.ipynb` or `seattle_ship.ipynb`

## Flooding between buildings experiment

After `make flow`, this creates a two-panel map of maximum depth and first flood
arrival on initially dry ground:

```sh
python street_flooding.py pike_place_surge12_fgout.nc
```

The script writes `street_flooding.png` and reports the number of
open-ground cells flooded in that view. Flooded>=0.2m, for cell that wasn't already wet at the start.

```sh
python street_flooding.py pike_place_surge12_fgout.nc --threshold 0.5 --view -122.356 -122.339 47.606 47.616
```
