from clawpack.geoclaw import fgout_tools
from clawpack.clawutil import data

def setrun(claw_pkg='geoclaw'):
    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)
    rundata = setgeo(rundata)
    clawdata = rundata.clawdata
    clawdata.num_dim = num_dim
    # domain: elliott bay (lon/lat)
    clawdata.lower[0] = -122.40
    clawdata.upper[0] = -122.32
    clawdata.lower[1] = 47.58
    clawdata.upper[1] = 47.64
    # 200x150 cells ~ 30m each
    clawdata.num_cells[0] = 200
    clawdata.num_cells[1] = 150
    # shallow water eqs: h, hu, hv
    clawdata.num_eqn = 3
    clawdata.num_aux = 3
    clawdata.capa_index = 2  # lon/lat capacity function
    clawdata.t0 = 0.0
    clawdata.restart = False

    # 20 full snapshots over 600s
    clawdata.output_style = 1
    clawdata.num_output_times = 20
    clawdata.tfinal = 600.0
    clawdata.output_t0 = True
    clawdata.output_format = 'ascii'
    clawdata.output_q_components = 'all'
    clawdata.output_aux_components = 'none'
    clawdata.output_aux_onlyonce = True
    clawdata.verbosity = 1

    # adaptive timestep via cfl condition
    clawdata.dt_variable = True
    clawdata.dt_initial = 0.1
    clawdata.dt_max = 1e+99
    clawdata.cfl_desired = 0.9
    clawdata.cfl_max = 1.0
    clawdata.steps_max = 50000

    # 2nd order shock-capturing scheme
    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2
    clawdata.num_waves = 3
    clawdata.limiter = ['mc', 'mc', 'mc']
    clawdata.use_fwaves = True
    clawdata.source_split = 'godunov'

    # waves leave the domain freely
    clawdata.num_ghost = 2
    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'
    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'
    clawdata.checkpt_style = 1

    # adaptive mesh refinement off
    amrdata = rundata.amrdata
    amrdata.amr_levels_max = 1  # higher vals turn on AMR
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
    # bunch of constants
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2  # lon/lat
    geo_data.earth_radius = 6367.5e3
    geo_data.coriolis_forcing = False
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3  # under 1mm treated as dry
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025  # open water
    geo_data.friction_depth = 20.0

    #AMR stuff
    refinement_data = rundata.refinement_data
    refinement_data.wave_tolerance = 1.e-1
    refinement_data.variable_dt_refinement_ratios = True

    # noaa dem from fetch_topo.py
    rundata.topo_data.topofiles.append([3, 'elliott_bay.tt3'])

    # 3m surface perturbation on west side drives the surge
    rundata.qinit_data.qinit_type = 4
    rundata.qinit_data.qinitfiles.append(['qinit.tt1'])

    # create the fgout grid for debris tracking later down the line: 300x225 every 5s
    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2
    fgout.output_format = 'binary32'
    fgout.nx = 300
    fgout.ny = 225
    fgout.x1 = -122.40
    fgout.x2 = -122.32
    fgout.y1 = 47.58
    fgout.y2 = 47.64
    fgout.tstart = 0.
    fgout.tend = 600.
    fgout.nout = 121
    rundata.fgout_data.fgout_grids.append(fgout)
    return rundata


if __name__ == '__main__':
    import sys
    rundata = setrun(*sys.argv[1:])
    rundata.write()
