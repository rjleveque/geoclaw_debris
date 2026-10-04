'''first attempt at ship tracking in geoclaw, seems to work at least in the small trial run notebook'''
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import MultiPoint, Polygon
from shapely.ops import transform

class Ship:
    def __init__(self, length_m=60., beam_m=12., draft_m=3., block_coefficient=.7, water_density=1025., drag_coefficient=1., under_keel_clearance_m=.5):
        self.length_m = length_m
        self.beam_m = beam_m
        self.draft_m = draft_m
        self.water_density = water_density
        self.drag_coefficient = drag_coefficient
        self.under_keel_clearance_m = under_keel_clearance_m
        self.mass_kg = water_density * length_m * beam_m * draft_m * block_coefficient
        self.yaw_inertia_kg_m2 = self.mass_kg * (length_m**2 + beam_m**2) / 12.

class LocalProjection:
    def __init__(self, lon0, lat0):
        self.lon0, self.lat0 = lon0, lat0
        phi = np.deg2rad(lat0)
        d = 1 - 0.00669437999014*np.sin(phi)**2
        self.m_lon = np.pi/180 * 6378137*np.cos(phi)/np.sqrt(d)
        self.m_lat = np.pi/180 * 6378137*(1-0.00669437999014)/d**1.5

    def to_xy(self, lon, lat): return (np.asarray(lon)-self.lon0)*self.m_lon, (np.asarray(lat)-self.lat0)*self.m_lat
    def to_lonlat(self, x, y): return np.asarray(x)/self.m_lon+self.lon0, np.asarray(y)/self.m_lat+self.lat0

class Flow(LocalProjection):
    def __init__(self, dataset, current_scale=1.):
        self.ds = dataset
        self.lon = np.asarray(dataset['lon'][:])
        self.lat = np.asarray(dataset['lat'][:])
        self.times = np.asarray(dataset['time'][:])
        self.current_scale = current_scale
        self.frames = {}
        super().__init__(self.lon.mean(), self.lat.mean())

    def frame(self, name, n):
        key = name, n
        if key not in self.frames: self.frames[key] = RegularGridInterpolator((self.lon, self.lat), np.asarray(self.ds[name][n]), bounds_error=False, fill_value=np.nan)
        return self.frames[key]

    def sample(self, name, xy, t):
        if t < self.times[0] or t > self.times[-1]: raise ValueError('time is outside the flow record')
        xy = np.atleast_2d(xy)
        lon, lat = self.to_lonlat(xy[:, 0], xy[:, 1])
        points = np.column_stack((lon, lat))
        n = min(np.searchsorted(self.times, t, side='right'), len(self.times)-1)
        m = max(n-1, 0)
        w = 0. if n == m else (t-self.times[m]) / (self.times[n]-self.times[m])
        values = (1-w)*self.frame(name, m)(points) + w*self.frame(name, n)(points)
        self.frames = {key: frame for key, frame in self.frames.items() if key[1] in (m, n)}
        return values*self.current_scale if name in ('u', 'v') else values

def hull_points(center, heading, ship):
    a, b = ship.length_m/2., ship.beam_m/2.
    xy = np.array([[0,0], [a,0], [-a,0], [a,b], [a,-b], [-a,b], [-a,-b], [0,b], [0,-b]])
    c, s = np.cos(heading), np.sin(heading)
    return np.asarray(center) + xy @ np.array([[c,s], [-s,c]])

def hull_polygon(center, heading, ship): return Polygon(hull_points(center, heading, ship)[[3,4,6,5]])

def check_pose(flow, buildings_utm, ship, center, heading, t):
    depths = flow.sample('h', hull_points(center, heading, ship), t)
    if not np.isfinite(depths).all(): return False, 'flow boundary', np.nan
    depth = float(depths.min())
    if depth < ship.draft_m + ship.under_keel_clearance_m: return False, 'grounding', depth
    if hull_polygon(center, heading, ship).intersects(buildings_utm): return False, 'building contact', depth
    return True, 'clear', depth

def simulate(flow, buildings_lonlat, ship, start_lonlat, heading_deg=0., start_time=0., end_time=600., dt=1.):
    if dt <= 0 or end_time <= start_time or not np.isclose((end_time-start_time)/dt, round((end_time-start_time)/dt)): raise ValueError('choose a positive dt that divides the time span')
    if min(ship.length_m, ship.beam_m, ship.draft_m, ship.mass_kg) <= 0: raise ValueError('ship dimensions and mass must be positive')

    buildings = transform(flow.to_xy, buildings_lonlat)
    center = np.array(flow.to_xy(*start_lonlat), dtype=float)
    heading = np.deg2rad(heading_deg)
    velocity = np.zeros(2)
    yaw_rate = 0.
    path = []
    stopped = None

    for t in np.arange(start_time, end_time+dt/2., dt):
        clear, reason, depth = check_pose(flow, buildings, ship, center, heading, t)
        if not clear and stopped is None: raise ValueError(f'invalid ship pose at t={t}: {reason}')
        path.append({'time_s': float(t), 'x_m': float(center[0]), 'y_m': float(center[1]), 'heading_deg': float(np.rad2deg(heading)), 'u_mps': float(velocity[0]), 'v_mps': float(velocity[1]), 'min_depth_m': depth, 'status': stopped or 'afloat'})
        if t >= end_time or stopped: continue

        points = hull_points(center, heading, ship)
        water = np.column_stack((flow.sample('u', points, t), flow.sample('v', points, t)))
        if not np.isfinite(water).all():
            stopped = 'flow boundary'
            velocity[:] = 0.
            yaw_rate = 0.
            continue

        offsets = points-center
        c, s = np.cos(heading), np.sin(heading)
        axes = np.array([[c,s], [-s,c]])
        point_velocity = velocity + yaw_rate*np.column_stack((-offsets[:,1], offsets[:,0]))
        relative = (water-point_velocity) @ axes.T
        area = ship.draft_m*np.array([ship.beam_m, ship.length_m]) / len(points)
        forces = (.5*ship.water_density*ship.drag_coefficient*area*relative*np.abs(relative)) @ axes
        force = forces.sum(axis=0)
        torque = np.sum(offsets[:,0]*forces[:,1] - offsets[:,1]*forces[:,0])
        next_velocity = velocity + dt*force/ship.mass_kg
        next_yaw_rate = yaw_rate + dt*torque/ship.yaw_inertia_kg_m2
        next_center = center + dt*(velocity+next_velocity)/2.
        next_heading = heading + dt*(yaw_rate+next_yaw_rate)/2.

        def try_step(fraction, center=center, heading=heading, t=t, points=points, next_center=next_center, next_heading=next_heading):
            position = center + fraction*(next_center-center)
            angle = heading + fraction*(next_heading-heading)
            valid, why, _ = check_pose(flow, buildings, ship, position, angle, t+fraction*dt)
            swept = MultiPoint(np.vstack((points, hull_points(position, angle, ship)))).convex_hull
            return valid and not swept.intersects(buildings), why, position, angle

        valid, reason, position, angle = try_step(1.)
        if valid:
            center, heading = position, angle
            velocity, yaw_rate = next_velocity, next_yaw_rate
            continue
        low, high = 0., 1.
        for _ in range(24): # binary search for how far the ship can go without collision/grounding
            fraction = (low+high)/2.
            valid, _, _, _ = try_step(fraction)
            if valid: low = fraction
            else: high = fraction
        _, _, center, heading = try_step(low)
        velocity[:] = 0.
        yaw_rate = 0.
        stopped = reason if reason != 'clear' else 'building contact'
    return path
