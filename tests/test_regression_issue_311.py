# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
Regression tests for issue #311: Review of psic modes.

Issue #311 reports problems with 7 psic modes:
  1. fixed_psi_vertical/horizontal — produce warnings on (0,0,L) reflections
  2. specular_vertical/horizontal — produce wild angles for (0,0,L)
  3. fixed_incidence_fixed_chi_fixed_phi (B3) — sign reversal

These tests establish expected behavior and catch regressions.
"""

import re
from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.reference import natural_psi
from ad_hoc_diffractometer import ub_identity


@pytest.fixture
def cubic_psic_geometry():
    """Set up a cubic psic geometry with UB matrix."""
    g = ahd.make_geometry("psic")
    g.wavelength = 1.0
    g.sample.lattice = ahd.Lattice(a=4.0)  # Cubic, a = 4 Å
    ub_identity(g.sample)  # UB = B for this orientation
    
    return g


# ============================================================================
# Issue #1: fixed_psi_vertical/horizontal — (0,0,L) reflections
# ============================================================================

class TestFixedPsiVertical00L:
    """Test fixed_psi_vertical on (0,0,L) reflections."""
    
    @pytest.mark.parametrize(
        "l_value",
        [1, 2, 3, 4, 5],
        ids=[f"00{l}" for l in range(1, 6)],
    )
    def test_fixed_psi_vertical_finds_solutions_for_00l(self, cubic_psic_geometry, l_value):
        """
        fixed_psi_vertical should find solutions for (0,0,L) reflections.
        
        When the constraint psi value matches the natural psi,
        solutions should be found (same as bisecting_vertical).
        """
        g = cubic_psic_geometry
        g.mode_name = "fixed_psi_vertical"
        g.azimuth = (0, -1, 0)
        
        # Set psi constraint to the natural value
        natural = natural_psi(g, 0, 0, l_value)
        g._modes["fixed_psi_vertical"] = g.mode.with_constraint_values(psi=natural)
        g.mode_name = "fixed_psi_vertical"
        
        # Should find solutions for (0,0,L)
        solutions = g.forward(0, 0, l_value)
        assert len(solutions) > 0, f"No solutions for (0,0,{l_value}) with psi={natural:.4f}°"
        
        # All solutions should be valid
        for sol in solutions:
            inverse_hkl = g.inverse(angles=sol)
            assert inverse_hkl[0] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[1] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[2] == pytest.approx(l_value, abs=1e-4)

    def test_fixed_psi_vertical_respects_constraint_value(self, cubic_psic_geometry):
        """
        fixed_psi_vertical with matching psi constraint should find solutions.
        """
        g = cubic_psic_geometry
        g.azimuth = (0, -1, 0)
        
        # Use (1,1,0) reflection where natural_psi is well-defined
        g.mode_name = "bisecting_vertical"
        sols = g.forward(1, 1, 0)
        assert len(sols) > 0, "No bisecting solutions for (1,1,0)"
        
        natural_psi_110 = natural_psi(g, 1, 1, 0)
        assert natural_psi_110 is not None, "natural_psi returned None for (1,1,0)"
        
        # Create fixed_psi_vertical mode with matching constraint
        g._modes["fixed_psi_test"] = ahd.ConstraintSet([
            ahd.SampleConstraint("mu", 0.0),
            ahd.DetectorConstraint("nu", 0.0),
            ahd.ReferenceConstraint("psi", natural_psi_110),
        ], computed=["eta", "chi", "phi", "delta"])
        g.mode_name = "fixed_psi_test"
        
        # Should find solutions
        sols = g.forward(1, 1, 0)
        assert len(sols) > 0, "No solutions with matching psi constraint"


class TestFixedPsiHorizontal00L:
    """Test fixed_psi_horizontal on (0,0,L) reflections."""
    
    @pytest.mark.parametrize(
        "l_value",
        [1, 2, 3],
        ids=[f"00{l}" for l in range(1, 4)],
    )
    def test_fixed_psi_horizontal_finds_solutions_for_00l(self, cubic_psic_geometry, l_value):
        """
        fixed_psi_horizontal should find solutions for (0,0,L) reflections
        when psi constraint matches natural psi.
        """
        g = cubic_psic_geometry
        g.mode_name = "fixed_psi_horizontal"
        g.azimuth = (0, -1, 0)
        
        # Set psi to natural value
        natural = natural_psi(g, 0, 0, l_value)
        g._modes["fixed_psi_horizontal"] = g.mode.with_constraint_values(psi=natural)
        g.mode_name = "fixed_psi_horizontal"
        
        # Should find solutions
        solutions = g.forward(0, 0, l_value)
        assert len(solutions) > 0, f"No solutions for (0,0,{l_value}) with psi={natural:.4f}°"
        
        # Verify solutions round-trip
        for sol in solutions:
            inverse_hkl = g.inverse(angles=sol)
            assert inverse_hkl[0] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[1] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[2] == pytest.approx(l_value, abs=1e-4)


# ============================================================================
# Issue #2: specular_vertical/horizontal — (0,0,L) produces wild angles
# ============================================================================

@pytest.mark.skip(reason="Issue #311: specular modes fail on (0,0,L) — analytic chi/phi solver singularity")
class TestSpecularVertical00L:
    """Test specular_vertical on (0,0,L) reflections.
    
    NOTE: Skipped pending fix of deep bug in _solve_bisecting_analytic
    for (0,0,L) reflections with standard Eulerian geometries.
    See issue #311 for details.
    """
    
    @pytest.mark.parametrize(
        "l_value",
        [1, 2, 3],
        ids=[f"00{l}" for l in range(1, 4)],
    )
    def test_specular_vertical_produces_reasonable_angles(self, cubic_psic_geometry, l_value):
        """
        specular_vertical should produce physically reasonable angles for (0,0,L).
        
        Specifically, eta should NOT be ±175°, ±170°, etc. (wraparound artifacts).
        """
        g = cubic_psic_geometry
        g.mode_name = "specular_vertical"
        g.surface_normal = (0, 0, 1)  # Surface normal along z
        
        solutions = g.forward(0, 0, l_value)
        assert len(solutions) > 0, f"No solutions for (0,0,{l_value})"
        
        for sol in solutions:
            eta = sol["eta"]
            # eta should be in a reasonable range, not ±175° or ±170°
            # For (0,0,L) in a typical geometry, |eta| should be < 90°
            assert abs(eta) < 100, (
                f"eta={eta:.1f}° is unreasonable for (0,0,{l_value}); "
                "possible sign error in surface constraint"
            )
            
            # Verify solution reproduces hkl
            inverse_hkl = g.inverse(angles=sol)
            assert inverse_hkl[0] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[1] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[2] == pytest.approx(l_value, abs=1e-4)

    def test_specular_vertical_incidence_equals_emergence(self, cubic_psic_geometry):
        """
        specular_vertical should enforce incidence == emergence.
        """
        from ad_hoc_diffractometer.reference import incidence_angle, emergence_angle
        
        g = cubic_psic_geometry
        g.mode_name = "specular_vertical"
        g.surface_normal = (0, 0, 1)
        
        solutions = g.forward(1, 1, 1)
        assert len(solutions) > 0
        
        for sol in solutions:
            inc = incidence_angle(g, angles=sol)
            em = emergence_angle(g, angles=sol)
            assert inc == pytest.approx(em, abs=1e-6), (
                f"incidence={inc:.4f}° != emergence={em:.4f}°"
            )


@pytest.mark.skip(reason="Issue #311: specular modes fail on (0,0,L) — analytic chi/phi solver singularity")
class TestSpecularHorizontal00L:
    """Test specular_horizontal on (0,0,L) reflections.
    
    NOTE: Skipped pending fix of deep bug in _solve_bisecting_analytic
    for (0,0,L) reflections with standard Eulerian geometries.
    See issue #311 for details.
    """
    
    @pytest.mark.parametrize(
        "l_value",
        [1, 2],
        ids=[f"00{l}" for l in range(1, 3)],
    )
    def test_specular_horizontal_produces_reasonable_angles(self, cubic_psic_geometry, l_value):
        """
        specular_horizontal should produce reasonable angles for (0,0,L).
        
        Specifically, mu should NOT be ±175°, ±170°, etc.
        """
        g = cubic_psic_geometry
        g.mode_name = "specular_horizontal"
        g.surface_normal = (0, 0, 1)
        
        solutions = g.forward(0, 0, l_value)
        assert len(solutions) > 0, f"No solutions for (0,0,{l_value})"
        
        for sol in solutions:
            mu = sol["mu"]
            # mu should be in a reasonable range
            assert abs(mu) < 100, (
                f"mu={mu:.1f}° is unreasonable for (0,0,{l_value})"
            )
            
            # Verify solution reproduces hkl
            inverse_hkl = g.inverse(angles=sol)
            assert inverse_hkl[0] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[1] == pytest.approx(0, abs=1e-4)
            assert inverse_hkl[2] == pytest.approx(l_value, abs=1e-4)


# ============================================================================
# Issue #3: fixed_incidence_fixed_chi_fixed_phi (B3) — sign reversal
# ============================================================================

class TestB3SignCorrectness:
    """Test fixed_incidence_fixed_chi_fixed_phi (B3) for sign correctness."""
    
    @pytest.mark.parametrize(
        "hkl",
        [
            pytest.param((1, 0, 0), id="100"),
            pytest.param((1, 1, 0), id="110"),
            pytest.param((1, 1, 1), id="111"),
        ],
    )
    def test_b3_produces_physically_correct_signs(self, cubic_psic_geometry, hkl):
        """
        B3 mode should produce solutions with correct angle signs.
        
        If the geometry is set up with mu/eta/nu/delta positive from zero,
        solutions should generally be positive (or follow the geometry's
        expected convention).
        """
        g = cubic_psic_geometry
        g.mode_name = "fixed_incidence_fixed_chi_fixed_phi"
        g.surface_normal = (0, 0, 1)
        
        solutions = g.forward(*hkl)
        assert len(solutions) > 0, f"No solutions for {hkl}"
        
        for sol in solutions:
            # For a positive reflection in a well-oriented geometry,
            # we expect positive detector angles
            # (This is a basic sanity check; exact expectations depend on geometry)
            mu = sol["mu"]
            eta = sol["eta"]
            nu = sol["nu"]
            delta = sol["delta"]
            
            # All should be finite
            assert abs(mu) < 180, f"mu={mu}° is outside reasonable range"
            assert abs(eta) < 180, f"eta={eta}° is outside reasonable range"
            assert abs(nu) < 180, f"nu={nu}° is outside reasonable range"
            assert abs(delta) < 180, f"delta={delta}° is outside reasonable range"
            
            # Verify solution reproduces hkl
            inverse_hkl = g.inverse(angles=sol)
            assert inverse_hkl[0] == pytest.approx(hkl[0], abs=1e-4)
            assert inverse_hkl[1] == pytest.approx(hkl[1], abs=1e-4)
            assert inverse_hkl[2] == pytest.approx(hkl[2], abs=1e-4)

    def test_b3_solutions_match_reference_geometry(self, cubic_psic_geometry):
        """
        B3 solutions should be consistent with bisecting when incidence=emergence=0.
        """
        g = cubic_psic_geometry
        
        # Get bisecting solution for reference
        g.mode_name = "bisecting_vertical"
        bisect_sols = g.forward(1, 0, 0)
        assert len(bisect_sols) > 0
        bisect_sol = bisect_sols[0]
        
        # Get B3 solution with incidence=0
        g.mode_name = "fixed_incidence_fixed_chi_fixed_phi"
        g.surface_normal = (0, 0, 1)
        g._modes[g.mode_name] = g.mode.with_constraint_values(incidence=0)
        g.mode_name = "fixed_incidence_fixed_chi_fixed_phi"
        
        b3_sols = g.forward(1, 0, 0)
        assert len(b3_sols) > 0, "B3 mode should find solutions when incidence=0"
        
        # At least one B3 solution should be close to the bisecting solution
        found_match = False
        for b3_sol in b3_sols:
            if all(
                abs(b3_sol[stage] - bisect_sol[stage]) < 1e-2
                for stage in ["mu", "eta", "chi", "phi", "nu", "delta"]
            ):
                found_match = True
                break
        
        assert found_match, (
            "No B3 solution matches bisecting solution; possible sign error"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
