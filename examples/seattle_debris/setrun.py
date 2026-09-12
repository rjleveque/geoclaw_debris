from clawpack.geoclaw import fgout_tools
from clawpack.clawutil import data
import numpy as np

# resolution control: define target cell size in meters, determines num_cells over the domain.
DX_METERS = 10.0

# domain for the simulation (currently elliott bay values) shared with make_qinit.py
X1, X2 = -122.40, -122.32
Y1, Y2 = 47.58, 47.64
#stack topo files in order of increasing resolution
TOPOFILES = [
    [3, 'seattle_19arcsec_mhw.tt3'], # 1/9" (~3m) full domain
    [3, '../../SeattleBuildings/SeattleBuildings1_19s.asc'], # 1/9" topo with building footprints stamped in (+20m height over footprints)
]

def setrun(claw_pkg='geoclaw'):
    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)
    rundata = setgeo(rundata)
    clawdata = rundata.clawdata
    clawdata.num_dim = num_dim

    clawdata.lower[0] = X1
    clawdata.upper[0] = X2
    clawdata.lower[1] = Y1
    clawdata.upper[1] = Y2

    # deriving num cells from DX_METERS
    lat_mid = 0.5 * (Y1 + Y2)
    m_per_deg_lon = 111000. * np.cos(np.radians(lat_mid))
    m_per_deg_lat = 111000.
    clawdata.num_cells[0] = max(10, int(round((X2-X1) * m_per_deg_lon / DX_METERS)))
    clawdata.num_cells[1] = max(10, int(round((Y2-Y1) * m_per_deg_lat / DX_METERS)))
    print(f'Grid: {clawdata.num_cells[0]} x {clawdata.num_cells[1]} cells '
          f'at ~{DX_METERS:.0f}m resolution')

    clawdata.num_eqn = 3
    clawdata.num_aux = 3
    clawdata.capa_index = 2
    clawdata.t0 = 0.0
    clawdata.restart = False

    clawdata.output_style = 1
    clawdata.num_output_times = 20
    clawdata.tfinal = 600.0
    clawdata.output_t0 = True
    clawdata.output_format = 'ascii'
    clawdata.output_q_components = 'all'
    clawdata.output_aux_components = 'none'
    clawdata.output_aux_onlyonce = True
    clawdata.verbosity = 1

    clawdata.dt_variable = True
    clawdata.dt_initial = 0.1
    clawdata.dt_max = 1e+99
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    clawdata.steps_max = 500000

    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2
    clawdata.num_waves = 3
    clawdata.limiter = ['mc', 'mc', 'mc']
    clawdata.use_fwaves = True
    clawdata.source_split = 'godunov'

    clawdata.num_ghost = 2
    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'
    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'
    clawdata.checkpt_style = 1

    amrdata = rundata.amrdata
    amrdata.amr_levels_max = 1
    amrdata.refinement_ratios_x = [4]
    amrdata.refinement_ratios_y = [4]
    amrdata.refinement_ratios_t = [4]
    amrdata.aux_type = ['center', 'capacity', 'yleft']
    amrdata.flag_richardson = False
    amrdata.flag2refine = True
    amrdata.regrid_interval = 3
    amrdata.regrid_buffer_width = 3
    amrdata.clustering_cutoff = 0.7
    amrdata.verbosity_regrid = 0

    return rundata


def setgeo(rundata):
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3
    geo_data.coriolis_forcing = False
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 20.0

    refinement_data = rundata.refinement_data
    refinement_data.wave_tolerance = 1.e-1
    refinement_data.variable_dt_refinement_ratios = True

    for topo_entry in TOPOFILES:
        rundata.topo_data.topofiles.append(topo_entry)

    rundata.qinit_data.qinit_type = 4
    rundata.qinit_data.qinitfiles.append(['qinit.tt1'])

    # fgout grid matches simulation domain and resolution
    lat_mid = 0.5 * (Y1 + Y2)
    m_per_deg_lon = 111000. * np.cos(np.radians(lat_mid))
    m_per_deg_lat = 111000.
    fgout_nx = max(10, int(round((X2-X1) * m_per_deg_lon / DX_METERS)))
    fgout_ny = max(10, int(round((Y2-Y1) * m_per_deg_lat / DX_METERS)))

    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2
    fgout.output_format = 'binary32'
    fgout.nx = fgout_nx
    fgout.ny = fgout_ny
    fgout.x1 = X1
    fgout.x2 = X2
    fgout.y1 = Y1
    fgout.y2 = Y2
    fgout.tstart = 0.
    fgout.tend = 600.
    fgout.nout = 121
    rundata.fgout_data.fgout_grids.append(fgout)
    return rundata


if __name__ == '__main__':
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
