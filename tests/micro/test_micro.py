"""Tests for the micro-turbojet extensions (sizing, build, metrics, clearance, cases, export)."""
import json
import os
import tempfile
import unittest
import warnings

import numpy as np

warnings.filterwarnings("ignore")

from microbladegen.micro import (
    AIR, Gas, size_impeller, build_impeller, wiesner_slip, glassman_blade_count,
    clearance_budget, clearance_effect, disk_rim_growth, thermal_growth, gap_adjusted_height,
    apply_tip_gap_geometry, free_vortex_angles, list_cases, build_case, export_case,
    wrap_to_cylinder, case_summary, impeller_flow_table,
)


class GasTest(unittest.TestCase):
    def test_air_properties(self):
        self.assertAlmostEqual(float(AIR.speed_of_sound(288.15)), 340.3, delta=0.3)
        self.assertAlmostEqual(float(AIR.viscosity(288.15)), 1.789e-5, delta=0.02e-5)

    def test_static_from_total_roundtrip(self):
        T, p, rho = AIR.static_from_total(300.0, 2.0e5, 0.0)
        self.assertAlmostEqual(float(T), 300.0)
        self.assertAlmostEqual(float(p), 2.0e5)
        self.assertAlmostEqual(float(rho), 2.0e5 / (AIR.R * 300.0))
        T, _, _ = AIR.static_from_total(300.0, 2.0e5, 200.0)
        self.assertAlmostEqual(float(T), 300.0 - 200.0**2 / (2 * AIR.cp))


class ImpellerSizingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = size_impeller(0.25, 110000, 3.2)

    def test_euler_work_closes(self):
        s = self.s
        # specific work = (1 + f_disk) * (sigma - phi tan(b2)) * U2^2, f_disk = 0.05 default
        mu = s.slip_factor - 0.25 * np.tan(np.radians(30.0))
        self.assertAlmostEqual(s.specific_work, 1.05 * mu * s.u2**2, delta=1.0)

    def test_work_meets_pressure_ratio(self):
        s = self.s
        dh_is = AIR.cp * s.Tt1 * (s.pressure_ratio ** ((AIR.gamma - 1) / AIR.gamma) - 1)
        self.assertAlmostEqual(s.specific_work * s.eta_tt, dh_is, delta=1.0)

    def test_tip_speed_radius_consistent(self):
        s = self.s
        self.assertAlmostEqual(s.u2, 2 * np.pi * s.rpm / 60 * s.r2, places=6)

    def test_tip_mach_definition(self):
        s = self.s
        self.assertAlmostEqual(s.tip_mach, s.u2 / float(AIR.speed_of_sound(s.Tt1)), places=9)

    def test_exit_continuity(self):
        s = self.s
        Tt2 = s.Tt1 + s.specific_work / AIR.cp
        c2 = np.hypot(s.cr2, s.ctheta2)
        _, _, rho2 = AIR.static_from_total(Tt2, s.Pt1 * s.pressure_ratio, c2)
        m = float(rho2) * s.cr2 * 2 * np.pi * s.r2 * s.b2 * (1 - 0.10)
        self.assertAlmostEqual(m, s.mass_flow, places=6)

    def test_inlet_continuity(self):
        s = self.s
        _, _, rho1 = AIR.static_from_total(s.Tt1, s.Pt1, s.cx1)
        area = np.pi * (s.r1s**2 - s.r1h**2) * (1 - 0.05)
        self.assertAlmostEqual(float(rho1) * s.cx1 * area, s.mass_flow, places=6)

    def test_inducer_flow_angle(self):
        s = self.s
        self.assertAlmostEqual(np.degrees(np.arctan(2 * np.pi * s.rpm / 60 * s.r1s / s.cx1)), s.beta1s_flow, places=6)

    def test_slip_and_blade_count(self):
        self.assertAlmostEqual(wiesner_slip(12, 0.0), 1 - 1 / 12**0.7)
        self.assertGreater(wiesner_slip(16, 30.0), wiesner_slip(8, 30.0))
        self.assertGreater(glassman_blade_count(70.0), glassman_blade_count(60.0) * 0.5)
        self.assertEqual(self.s.n_total % 2, 0)

    def test_scaling_trends(self):
        hi = size_impeller(0.25, 130000, 3.2)
        self.assertLess(hi.r2, self.s.r2)               # same work, more rpm -> smaller wheel
        big = size_impeller(0.35, 110000, 3.2)
        self.assertGreater(big.r1s, self.s.r1s)         # more flow -> bigger inducer
        pr = size_impeller(0.25, 110000, 4.0)
        self.assertGreater(pr.u2, self.s.u2)            # more pressure ratio -> more tip speed

    def test_micro_values_realistic(self):
        s = self.s
        self.assertTrue(1.1 < s.tip_mach < 1.5)
        self.assertTrue(0.8 < s.inlet_rel_mach_shroud < 1.1)
        self.assertTrue(1.4 < s.diffusion_ratio < 1.95)
        self.assertTrue(3e-3 < s.b2 < 8e-3)
        self.assertTrue(5e5 < s.machine_reynolds < 3e6)

    def test_warning_notes(self):
        s = size_impeller(0.15, 150000, 3.0, backsweep=25.0, beta1s_flow=57.0)
        self.assertTrue(any("Diffusion ratio" in n for n in s.notes))


class ImpellerBuildTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = size_impeller(0.25, 110000, 3.2)
        cls.imp = build_impeller(cls.s, tip_clearance=0.3e-3)

    def test_angles_hit_targets(self):
        m = self.imp.measured
        for got, want in zip(m["inlet_angle_hub_mid_tip"], m["inlet_angle_target_hub_mid_tip"]):
            self.assertAlmostEqual(got, want, delta=1.5)
        for got in m["exit_angle_hub_mid_tip"]:
            self.assertAlmostEqual(got, m["exit_angle_target"], delta=1.5)

    def test_exit_radius(self):
        self.assertAlmostEqual(self.imp.measured["exit_radius"], self.s.r2, delta=2e-4)

    def test_blade_leans_one_way(self):
        """No S-curve: blade angle keeps one sign over the chord at every span."""
        from microbladegen.micro.impeller import _mean_camber_angles
        mb = self.imp.centrif.mainblade
        for j in (0, mb.ss_cyl_pts.shape[0] // 2, mb.ss_cyl_pts.shape[0] - 1):
            ang, _ = _mean_camber_angles(mb, j)
            self.assertTrue(np.all(ang[2:-2] < 0))

    def test_thickness_close_to_request(self):
        t = self.imp.measured["max_thickness_hub_mid_tip"]
        self.assertAlmostEqual(t[-1], 0.7e-3, delta=0.2e-3)
        self.assertGreater(t[0], t[-1])

    def test_blade_counts(self):
        self.assertEqual(self.imp.nblades, 6)
        self.assertIsNotNone(self.imp.centrif.splitterblade)

    def test_tip_gap_geometry(self):
        gap = 0.3e-3
        mb = self.imp.centrif.mainblade
        tip = mb.ss_cyl_pts[-1]
        sh = self.imp.shroud
        dense = np.linspace(0, 1, 4000)
        xs = np.interp(dense, np.linspace(0, 1, len(sh)), sh[:, 0])
        rs = np.interp(dense, np.linspace(0, 1, len(sh)), sh[:, 1])
        d = np.array([np.min(np.hypot(xs - x, rs - r)) for x, r in zip(tip[5:-5, 0], tip[5:-5, 1])])
        self.assertLessEqual(d.max(), gap * 1.02)
        self.assertGreater(d.min(), 0.4 * gap)

    def test_cart_matches_cyl(self):
        mb = self.imp.centrif.mainblade
        r, th = mb.ss_cyl_pts[:, :, 1], mb.ss_cyl_pts[:, :, 2]
        np.testing.assert_allclose(mb.ss_cart_pts[:, :, 1], r * np.sin(th), atol=1e-12)
        np.testing.assert_allclose(mb.ss_cart_pts[:, :, 2], r * np.cos(th), atol=1e-12)

    def test_zero_gap_leaves_tip_on_shroud(self):
        imp0 = build_impeller(self.s, tip_clearance=0.0, refine_passes=0)
        tip = imp0.centrif.mainblade.ss_cyl_pts[-1]
        sh = imp0.shroud
        dense = np.linspace(0, 1, 4000)
        xs = np.interp(dense, np.linspace(0, 1, len(sh)), sh[:, 0])
        rs = np.interp(dense, np.linspace(0, 1, len(sh)), sh[:, 1])
        d = np.array([np.min(np.hypot(xs - x, rs - r)) for x, r in zip(tip[5:-5, 0], tip[5:-5, 1])])
        self.assertLess(d.max(), 2e-5)


class TipClearanceTest(unittest.TestCase):
    def test_disk_growth_closed_form(self):
        r, u = 0.038, 440.0
        m = dict(rho=4430.0, E=113.8e9, nu=0.342)
        expected = (1 - m["nu"]) * m["rho"] * u**2 * r / (4 * m["E"])
        self.assertAlmostEqual(disk_rim_growth(r, u, "titanium_6al4v"), expected, places=12)

    def test_growth_scales(self):
        self.assertAlmostEqual(disk_rim_growth(0.04, 400, "aluminium_7075") / disk_rim_growth(0.04, 200, "aluminium_7075"), 4.0)
        self.assertAlmostEqual(thermal_growth(0.04, 100, "inconel_718"), 13.0e-6 * 100 * 0.04)

    def test_budget_ordering(self):
        b = clearance_budget(0.038, 440.0)
        self.assertLessEqual(b.rss, b.worst_case)
        self.assertGreaterEqual(b.recommended, b.disk_growth + b.thermal_mismatch)
        self.assertGreater(b.recommended, 0.05e-3)
        hot = clearance_budget(0.038, 440.0, rotor_delta_T=300.0, casing_delta_T=20.0)
        self.assertGreater(hot.thermal_mismatch, b.thermal_mismatch)

    def test_effect_linear(self):
        a = clearance_effect(0.1e-3, 5e-3, "centrifugal")
        b = clearance_effect(0.2e-3, 5e-3, "centrifugal")
        self.assertAlmostEqual(b.delta_efficiency_points, 2 * a.delta_efficiency_points)
        self.assertAlmostEqual(a.ratio_percent, 2.0)
        self.assertEqual(clearance_effect(0.05e-3, 10e-3, "axial_compressor").status, "ok")
        self.assertEqual(clearance_effect(0.5e-3, 5e-3, "axial_compressor").status, "high")

    def test_gap_adjusted_height(self):
        self.assertAlmostEqual(gap_adjusted_height(0.014, 0.0002), 0.0138)
        with self.assertRaises(ValueError):
            gap_adjusted_height(0.001, 0.002)


class AxialMetricsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rotor = build_case("axial_rotor_064krpm")
        cls.ngv = build_case("ngv_micro_turbojet")

    def test_free_vortex_rotor_angles(self):
        a = self.rotor.angles
        self.assertAlmostEqual(a["beta1"].iloc[1], 52.83, delta=0.1)
        self.assertAlmostEqual(a["beta2"].iloc[1], 37.98, delta=0.1)
        self.assertLess(a["beta2"].iloc[0], a["beta1"].iloc[0])      # compressor rotor diffuses

    def test_stator_angles_equal_absolute(self):
        a = self.ngv.angles
        np.testing.assert_allclose(a["beta2"], a["alpha2"])

    def test_solidity_and_pitch(self):
        m = self.rotor.metrics
        np.testing.assert_allclose(m["solidity"], m["chord_mm"] / m["pitch_mm"])
        np.testing.assert_allclose(m["pitch_mm"], 2 * np.pi * m["radius_mm"] / 29)

    def test_tip_blade_mach(self):
        m = self.rotor.metrics.iloc[-1]
        self.assertGreater(m["tip_blade_mach"], 0.9)
        self.assertLess(m["tip_blade_mach"], 1.2)
        self.assertGreater(m["tip_blade_mach"], m["inlet_mach"])     # pre-swirl lowers relative velocity

    def test_de_haller_only_for_diffusing_rows(self):
        self.assertTrue(np.all(self.rotor.metrics["de_haller"] > 0.7))
        self.assertTrue(self.ngv.metrics["de_haller"].isna().all())
        self.assertTrue((self.ngv.metrics["velocity_ratio"] > 1.0).all())

    def test_ngv_exit_mach_subsonic_and_hub_highest(self):
        m = self.ngv.metrics["exit_mach"]
        self.assertTrue((m < 1.0).all())
        self.assertGreater(m.iloc[0], m.iloc[-1])

    def test_hot_gas_used(self):
        self.assertEqual(self.ngv.flow.gas.gamma, 1.33)

    def test_reynolds_positive(self):
        self.assertTrue((self.rotor.metrics["reynolds_chord"] > 5e4).all())
        self.assertTrue((self.ngv.metrics["reynolds_chord_exit"] > 5e4).all())

    def test_blade_height_reflects_gap(self):
        zz = self.rotor.blade.zz
        self.assertAlmostEqual(zz[-1], (52.3 - 38.7 - 0.2) * 1e-3, places=9)


class WavyTest(unittest.TestCase):
    """The wavy stretch acts on the suction side and scales linearly with the ratio (documented behaviour)."""

    @classmethod
    def setUpClass(cls):
        from microbladegen.aero import AirfoilWavy
        from microbladegen.helper import StackType
        from microbladegen.micro import load_case, axial_row_foils
        foils, spans, _, height, *_ = axial_row_foils(load_case("ngv_micro_turbojet"))
        cls.n = 30
        cls.t = np.linspace(0, 2 * np.pi * 3, cls.n)

        def build(amp):
            b = AirfoilWavy(profileArray=foils, profile_loc=list(spans), height=height)
            b.stack(StackType.Centroid)
            b.build(nProfiles=cls.n, num_points=80, trailing_edge_points=20)
            if amp:
                z = np.zeros(cls.n)
                b.stretch_thickness_chord_te(amp * np.sin(cls.t), z, 0.5 * amp * np.sin(cls.t), z, z, z)
            return b
        cls.base, cls.w1, cls.w2 = build(0.0), build(0.05), build(0.10)

    def _disp(self, w):
        return np.array([np.max(np.linalg.norm(w.shft_ss[j] - self.base.shft_ss[j], axis=1)) for j in range(self.n)])

    def test_displacement_linear_in_ratio(self):
        d1, d2 = self._disp(self.w1), self._disp(self.w2)
        self.assertGreater(d1.max(), 5e-5)
        self.assertAlmostEqual(d2.max() / d1.max(), 2.0, delta=0.1)

    def test_displacement_oscillates_along_span(self):
        d = self._disp(self.w2)
        interior = d[1:-1]
        peaks = np.sum((interior[1:-1] > interior[:-2]) & (interior[1:-1] > interior[2:]))
        self.assertGreaterEqual(peaks, 2)


class CasesTest(unittest.TestCase):
    def test_list(self):
        names = list_cases()
        for n in ("impeller_250g", "impeller_150g", "axial_rotor_064krpm", "ngv_micro_turbojet"):
            self.assertIn(n, names)

    def test_case_files_are_valid_json(self):
        from importlib import resources
        folder = resources.files("microbladegen.micro").joinpath("cases")
        for p in folder.iterdir():
            if p.name.endswith(".json"):
                d = json.loads(p.read_text())
                self.assertIn("type", d)
                self.assertIn("description", d)

    def test_unknown_case(self):
        with self.assertRaises(FileNotFoundError):
            build_case("nope")

    def test_impeller_cases_build(self):
        for n in ("impeller_250g", "impeller_150g"):
            c = build_case(n)
            self.assertEqual(len(c.metrics), 3)
            self.assertGreater(c.clearance_budget.recommended, 0)
            self.assertLessEqual(c.sizing.diffusion_ratio, 1.9)

    def test_summary_text(self):
        txt = case_summary(build_case("impeller_250g"))
        for key in ("MEAN-LINE", "SECTIONS", "TIP CLEARANCE", "Tip Mach"):
            self.assertIn(key, txt)
        txt = case_summary(build_case("ngv_micro_turbojet"))
        for key in ("FLOW ANGLES", "SECTION METRICS", "NOTES"):
            self.assertIn(key, txt)


class ExportTest(unittest.TestCase):
    def test_wrap_to_cylinder(self):
        p = np.array([[0.01, 0.0, 0.0], [0.01, 0.02, 0.01]])
        w = wrap_to_cylinder(p, 0.03)
        self.assertAlmostEqual(float(np.hypot(w[0, 1], w[0, 2])), 0.03)
        self.assertAlmostEqual(float(np.hypot(w[1, 1], w[1, 2])), 0.04)
        self.assertAlmostEqual(float(np.arctan2(w[1, 1], w[1, 2])), 0.02 / 0.04)

    def test_axial_export(self):
        case = build_case("ngv_micro_turbojet")
        with tempfile.TemporaryDirectory() as d:
            man = export_case(case, d)
            self.assertAlmostEqual(man["periodic_sector_deg"], 360 / 17)
            for f in ("manifest.json", "blade_surface.stl", "summary_sections.csv", "summary_angles.csv"):
                self.assertTrue(os.path.exists(os.path.join(d, f)), f)
            first = np.loadtxt(os.path.join(d, "sections", man["section_files"]["blade"][0]))
            r = np.hypot(first[:, 1], first[:, 2])
            self.assertAlmostEqual(r.min(), 0.028, delta=2e-5)
            last = np.loadtxt(os.path.join(d, "sections", man["section_files"]["blade"][-1]))
            self.assertAlmostEqual(np.hypot(last[:, 1], last[:, 2]).max(), 0.038, delta=2e-5)
            from stl import mesh
            m = mesh.Mesh.from_file(os.path.join(d, "blade_surface.stl"))
            self.assertGreater(len(m.vectors), 1000)
            rr = np.hypot(m.vectors[..., 1], m.vectors[..., 2])
            self.assertLessEqual(rr.max(), 0.038 + 1e-5)
            self.assertGreaterEqual(rr.min(), 0.028 - 1e-5)
            with open(os.path.join(d, "manifest.json")) as f:
                self.assertEqual(json.load(f)["n_blades"], 17)

    def test_impeller_export(self):
        case = build_case("impeller_250g")
        with tempfile.TemporaryDirectory() as d:
            man = export_case(case, d)
            self.assertAlmostEqual(man["periodic_sector_deg"], 60.0)
            self.assertIn("splitter", man["section_files"])
            for f in ("hub_curve.csv", "shroud_curve.csv", "main_blade_surface.stl", "splitter_blade_surface.stl",
                      "summary_flow.csv", "summary_sections.csv", "manifest.json"):
                self.assertTrue(os.path.exists(os.path.join(d, f)), f)
            tip = np.loadtxt(os.path.join(d, "sections", man["section_files"]["main"][-1]))
            r = np.hypot(tip[:, 1], tip[:, 2])
            self.assertAlmostEqual(r.max(), case.sizing.r2, delta=2e-4)
            self.assertAlmostEqual(r.min(), case.sizing.r1s, delta=5e-4)
            self.assertAlmostEqual(man["tip_clearance_m"], case.clearance_budget.recommended)


if __name__ == "__main__":
    unittest.main()
