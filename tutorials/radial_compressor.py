"""Radial (Centripetal) Compressor Example

Builds a centripetal compressor geometry with a splitter blade using MicroBladeGen.
Demonstrates endwall construction, blade + splitter building, plotting, and data export.
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
    x_stretch_factor = 1.3
    rhub_out = 0.009  # meters
    nblades = 8

    # --- Build Endwalls (compressor: 270 to 360 degrees) ---
    hub1, shroud1 = build_endwalls(
        radius=radius,
        inlet_hub_shroud_ratio=inlet_hub_shroud_ratio,
        outlet_hub_shroud_ratio=outlet_hub_shroud_ratio,
        x_stretch_factor=x_stretch_factor,
        rhub_out=rhub_out,
        alpha_start=270, alpha_stop=360,
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

    # --- Build Centripetal Compressor ---
    cen = Centrif(blade_position=(0.0, 1.0), use_mid_wrap_angle=True, use_ray_camber=False)
    cen.add_hub(hub2[:, 0], hub2[:, 1])
    cen.add_shroud(shroud2[:, 0], shroud2[:, 1])

    TE_Cut = False
    te_props = TrailingEdgeProperties(TE_Cut=TE_Cut, TE_Radius=0.05)

    hub = CentrifProfile(percent_span=0, LE_Thickness=0.05,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=-30,
                         LE_Metal_Angle_Loc=0.1, TE_Metal_Angle_Loc=0.9,
                         ss_thickness=[0.05, 0.05, 0.05, 0.05, 0.05],
                         ps_thickness=[0.05, 0.05, 0.05, 0.05, 0.05],
                         wrap_angle=-25, use_ray_camber=True)
                        #  wrap_displacements=[0.5, 0.5, 0.5, 0.5],
                        #  wrap_displacement_locs=[0.3, 0.4, 0.5, 0.7])

    te_props = TrailingEdgeProperties(TE_Cut=TE_Cut, TE_Radius=0.05)

    mid = CentrifProfile(percent_span=0.5, LE_Thickness=0.03,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=-30,
                         LE_Metal_Angle_Loc=0.1, TE_Metal_Angle_Loc=0.9,
                         ss_thickness=[0.05, 0.05, 0.05, 0.05, 0.05],
                         ps_thickness=[0.05, 0.05, 0.05, 0.05, 0.05],
                         wrap_angle=-25, use_ray_camber=True)
                        #  wrap_displacements=[0.5, 0.5, 0.5, 0.5],
                        #  wrap_displacement_locs=[0.3, 0.4, 0.5, 0.7])

    tip = CentrifProfile(percent_span=1, LE_Thickness=0.03,
                         trailing_edge_properties=te_props,
                         LE_Metal_Angle=0, TE_Metal_Angle=-30,
                         LE_Metal_Angle_Loc=0.1, TE_Metal_Angle_Loc=0.9,
                         ss_thickness=[0.05, 0.05, 0.05, 0.05],
                         ps_thickness=[0.05, 0.05, 0.05, 0.05],
                         wrap_angle=-25, use_ray_camber=True)
                        #  wrap_displacements=[0.5, 0.5, 0.5, 0.5],
                        #  wrap_displacement_locs=[0.3, 0.4, 0.5, 0.7])

    cen.add_profile(hub)
    cen.add_profile(mid)
    cen.add_profile(tip)
    cen.add_splitter([hub, mid, tip], splitter_start=0.5)

    cen.build(npts_span=10, npts_chord=100, nblades=nblades)

    # --- Plots ---
    # cen.plot_camber()
    # cen.plot_mp_profile()
    # cen.plot()
    cen.plot_fullwheel()

    # --- Export Data ---
    # Mainblade
    write_data('centrif_compressor_mainblade_xyz_pts.csv',
               cen.mainblade.ss_cart_pts, cen.mainblade.ps_cart_pts,
               'ss_x,ss_y,ss_z,ps_x,ps_y,ps_z')

    write_data('centrif_compressor_mainblade_xrth.csv',
               cen.mainblade.ss_cyl_pts, cen.mainblade.ps_cyl_pts,
               'ss_x,ss_r,ss_th,ps_x,ps_r,ps_th')

    write_data_mp('centrif_compressor_mainblade_mp.csv',
                  cen.mainblade.ss_mp_pts, cen.mainblade.ps_mp_pts,
                  'ss_mp,ss_theta,ps_mp,ps_theta')

    # Splitter
    write_data('centrif_compressor_splitter_xyz_pts.csv',
               cen.splitterblade.ss_cart_pts, cen.splitterblade.ps_cart_pts,
               'ss_x,ss_y,ss_z,ps_x,ps_y,ps_z')

    write_data('centrif_compressor_splitter_xrth.csv',
               cen.splitterblade.ss_cyl_pts, cen.splitterblade.ps_cyl_pts,
               'ss_x,ss_r,ss_th,ps_x,ps_r,ps_th')

    write_data_mp('centrif_compressor_splitter_mp.csv',
                  cen.splitterblade.ss_mp_pts, cen.splitterblade.ps_mp_pts,
                  'ss_mp,ss_theta,ps_mp,ps_theta')
