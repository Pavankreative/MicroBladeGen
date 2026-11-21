import numpy as np
import pytest
from microbladegen.helper.pspline import pspline
from microbladegen.helper.resample import resample_curve, resample_by_curvature, arc_length_param_strict


def test_spline_preserves_traversal_and_query_order():
    spline = pspline([3, 2, 2, 0], [0, 0, 0, 0])
    points, _ = spline([1., 0., .5])
    np.testing.assert_allclose(points, [[0, 0], [3, 0], [1.5, 0]], atol=1e-6)


@pytest.mark.parametrize('function', [resample_curve, resample_by_curvature])
def test_closed_curve_keeps_closing_endpoint(function):
    points = np.array([[0, 0], [1, 0], [1, 0], [1, 1], [0, 0]])
    out = function(points, 21)
    np.testing.assert_allclose(out[[0, -1]], [[0, 0], [0, 0]], atol=1e-10)


@pytest.mark.parametrize('points', [[[2, 0], [0, 0]], [[2, 0, 0], [1, 0, 0], [0, 0, 0]]])
def test_short_curves(points):
    out = resample_curve(points, 9)
    np.testing.assert_allclose(out[[0, -1]], np.array(points)[[0, -1]], atol=1e-12)


@pytest.mark.parametrize('points', [[[0, 0], [0, 0]], [[0, 0], [np.nan, 1]], [[0, 0]]])
def test_invalid_geometry(points):
    with pytest.raises(ValueError):
        arc_length_param_strict(points)


@pytest.mark.parametrize('count', [0, 1, 2.5, True])
def test_invalid_count(count):
    with pytest.raises(ValueError):
        resample_curve([[0, 0], [1, 0]], count)


@pytest.mark.parametrize('z', [None, [0, 0, 0]])
def test_three_point_cubic_spline(z):
    from microbladegen.helper.pspline import spline_type
    spline = pspline([2, 1, 0], [0, 0, 0], z, method=spline_type.spline)
    points, _ = spline([0, .5, 1])
    expected = np.array([[2, 0], [1, 0], [0, 0]])
    if z is not None:
        expected = np.column_stack([expected, np.zeros(3)])
    np.testing.assert_allclose(points, expected, atol=1e-9)
