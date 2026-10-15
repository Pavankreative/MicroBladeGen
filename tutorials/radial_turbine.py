"""Radial Turbine Example

Builds a radial turbine geometry using MicroBladeGen with hub, mid, and tip profiles.
Demonstrates endwall construction, blade building, plotting, and data export.
"""
from microbladegen.aero import Centrif, CentrifProfile, TrailingEdgeProperties
from microbladegen.helper import arc
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt


def build_endwalls(radius: float, inlet_hub_shroud_ratio: float,
                   outlet_hub_shroud_ratio: float, x_stretch_factor: float,
                   rhub_out: float, alpha_start: float = 180, alpha_stop: float = 270):
    shroud = arc(xc=0, yc=0, radius=radius, alpha_start=alpha_start, alpha_stop=alpha_stop)
    hub = arc(xc=0, yc=0, radius=radius / inlet_hub_shroud_ratio, alpha_start=alpha_start, alpha_stop=alpha_stop)

    [xhub, rhub] = hub.get_point(np.linspace(0, 1, 100))
    [xshroud, rshroud] = shroud.get_point(np.linspace(0, 1, 100))

    rhub = rhub / outlet_hub_shroud_ratio
    xhub *= x_stretch_factor
    xshroud *= x_stretch_factor

    hub = np.vstack([xhub, rhub]).transpose()
    shroud = np.vstack([xshroud, rshroud]).transpose()
    shroud[:, 1] += -hub[:, 1].min() + rhub_out
    hub[:, 1] += -hub[:, 1].min() + rhub_out
    return hub, shroud


def compute_normals(x: npt.NDArray, y: npt.NDArray):
    dx = np.gradient(x)
    dy = np.gradient(y)
    length = np.hypot(dx, dy)
    nx = -dy / length
    ny = dx / length
    return nx, ny


def offset_curve(x: npt.NDArray, y: npt.NDArray, offset_distance: float):
    nx, ny = compute_normals(x, y)
    x_offset = x + offset_distance * nx
    y_offset = y + offset_distance * ny
    return x_offset, y_offset


def write_data(filename: str, ss_pts: npt.NDArray, ps_pts: npt.NDArray, headers: str):
    with open(filename, 'w') as f:
        f.write(f'span_index,{headers}\n')
        for i in range(ss_pts.shape[0]):
            for j in range(ss_pts.shape[1]):
                x1, y1, z1 = ss_pts[i, j]
                x2, y2, z2 = ps_pts[i, j]
                f.write(f"{i:d}, {x1:.8f}, {y1:.8f}, {z1:.8f}, {x2:.8f}, {y2:.8f}, {z2:.8f}\n")


def write_data_mp(filename: str, ss_pts: npt.NDArray, ps_pts: npt.NDArray, headers: str):
    with open(filename, 'w') as f:
        f.write(f'span_index,{headers}\n')
        for i in range(ss_pts.shape[0]):
            for j in range(ss_pts.shape[1]):
                x1, y1 = ss_pts[i, j]
                x2, y2 = ps_pts[i, j]
                f.write(f"{i:d}, {x1:.8f}, {y1:.8f}, {x2:.8f}, {y2:.8f}\n")


if __name__ == '__main__':
    # --- Endwall Parameters ---
    radius = 0.04  # meters
    inlet_hub_shroud_ratio = 0.85
    outlet_hub_shroud_ratio = 0.8
    x_stretch_factor = 1.1
    rhub_out = 0.009  # meters
    nblades = 8

    # --- Build Endwalls ---
    hub1, shroud1 = build_endwalls(
        radius=radius,
        inlet_hub_shroud_ratio=inlet_hub_shroud_ratio,
        outlet_hub_shroud_ratio=outlet_hub_shroud_ratio,
        x_stretch_factor=x_stretch_factor,
        rhub_out=rhub_out,
    )
    hub2 = np.vstack(offset_curve(hub1[:, 0], hub1[:, 1], -radius * 0.02)).T
    shroud2 = np.vstack(offset_curve(shroud1[:, 0], shroud1[:, 1], radius * 0.02)).T

    # Plot endwalls
    plt.figure(num=1, clear=True, dpi=150)
    plt.plot(hub1[:, 0], hub1[:, 1], 'k', linewidth=1.5, label='hub-1')
    plt.plot(shroud1[:, 0], shroud1[:, 1], 'k', linewidth=1.5, label='shroud-1')
    plt.plot(hub2[:, 0], hub2[:, 1], 'm', linewidth=1.5, label='hub-2')
    plt.plot(shroud2[:, 0], shroud2[:, 1], 'm', linewidth=1.5, label='shroud-2')
    plt.xlabel('x-axial')
    plt.ylabel('radius')
    plt.axis('equal')
    plt.legend()
    plt.show()

    # --- Build Radial Turbine ---
    rturbine = Centrif(blade_position=(0.0, 1.0), use_mid_wrap_angle=True, use_ray_camber=True)
    rturbine.add_hub(hub2[:, 0], hub2[:, 1])
    rturbine.add_shroud(shroud2[:, 0], shroud2[:, 1])

    TE_Cut = False
    te_props = TrailingEdgeProperties(TE_Cut=TE_Cut, TE_Radius=0.06)

    hub = CentrifProfile(percent_span=0, LE_Thickness=0.10,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=60,
                         LE_Metal_Angle_Loc=0.2, TE_Metal_Angle_Loc=0.8,
                         ss_thickness=[0.12, 0.12, 0.12, 0.12],
                         ps_thickness=[0.12, 0.12, 0.12, 0.12],
                         wrap_angle=-20,
                         wrap_displacements=[0, 0],
                         wrap_displacement_locs=[0, 0])

    mid = CentrifProfile(percent_span=0.5, LE_Thickness=0.06,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=60,
                         LE_Metal_Angle_Loc=0.2, TE_Metal_Angle_Loc=0.8,
                         ss_thickness=[0.07, 0.07, 0.07, 0.05],
                         ps_thickness=[0.07, 0.07, 0.07, 0.05],
                         wrap_angle=-20,
                         wrap_displacements=[0, 0],
                         wrap_displacement_locs=[0, 0])

    tip = CentrifProfile(percent_span=1, LE_Thickness=0.06,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=60,
                         LE_Metal_Angle_Loc=0.2, TE_Metal_Angle_Loc=0.8,
                         ss_thickness=[0.05, 0.05, 0.05, 0.05],
                         ps_thickness=[0.05, 0.05, 0.05, 0.05],
                         wrap_angle=-20,
                         wrap_displacements=[0, 0],
                         wrap_displacement_locs=[0, 0])

    rturbine.add_profile(hub)
    rturbine.add_profile(mid)
    rturbine.add_profile(tip)

    rturbine.build(npts_span=20, npts_chord=100, nblades=nblades)

    # --- Plots ---
    rturbine.plot_camber()
    rturbine.plot_mp_profile()
    rturbine.plot()
    rturbine.plot_fullwheel()

    # --- Export Data ---
    write_data('rturbine_xyz_pts.csv',
               rturbine.mainblade.ss_cart_pts, rturbine.mainblade.ps_cart_pts,
               'ss_x,ss_y,ss_z,ps_x,ps_y,ps_z')

    write_data('rturbine_xrth.csv',
               rturbine.mainblade.ss_cyl_pts, rturbine.mainblade.ps_cyl_pts,
               'ss_x,ss_r,ss_th,ps_x,ps_r,ps_th')

    write_data_mp('rturbine_mp.csv',
                  rturbine.mainblade.ss_mp_pts, rturbine.mainblade.ps_mp_pts,
                  'ss_mp,ss_theta,ps_mp,ps_theta')
