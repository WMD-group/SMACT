from __future__ import annotations

import pytest

from smact.lattice_parameters import (
    b2,
    b10,
    bcc,
    bct,
    cubic_perovskite,
    diamond,
    fcc,
    hcp,
    rocksalt,
    stuffed_wurtzite,
    wurtzite,
    zincblende,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

LatticeResult = tuple[float, float, float, float, float, float]
CUBIC_ANGLES = (90, 90, 90)
HEX_ANGLES = (90, 90, 120)


def _assert_positive_lengths(result: LatticeResult):
    a, b, c, _alpha, _beta, _gamma = result
    assert a > 0
    assert b > 0
    assert c > 0


def _assert_six_tuple(result: LatticeResult):
    assert isinstance(result, tuple)
    assert len(result) == 6
    for v in result:
        assert isinstance(v, (int, float))


def _assert_angles(result: LatticeResult, expected: tuple[int, int, int]):
    _, _, _, alpha, beta, gamma = result
    assert alpha == expected[0]
    assert beta == expected[1]
    assert gamma == expected[2]


# ---------------------------------------------------------------------------
# Tests for functions taking a list of shannon radii
# ---------------------------------------------------------------------------


class TestCubicPerovskite:
    @pytest.mark.parametrize("radii", [[1.0, 0.6, 1.4], [1.5, 0.5, 1.2], [0.8, 0.4, 0.9]])
    def test_returns_valid_tuple(self, radii):
        result = cubic_perovskite(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = cubic_perovskite([1.0, 0.6, 1.4])
        a, b, c, *_ = result
        assert a == b == c


class TestWurtzite:
    @pytest.mark.parametrize("radii", [[1.0, 0.5], [0.7, 0.7], [1.5, 0.3]])
    def test_returns_valid_tuple(self, radii):
        result = wurtzite(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, HEX_ANGLES)

    def test_hexagonal_symmetry(self):
        result = wurtzite([1.0, 0.5])
        a, b, _c, *_ = result
        assert a == b


class TestRocksalt:
    @pytest.mark.parametrize("radii", [[1.0, 1.0], [0.5, 1.5], [1.2, 0.8]])
    def test_returns_valid_tuple(self, radii):
        result = rocksalt(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = rocksalt([1.0, 1.0])
        a, b, c, *_ = result
        assert a == b == c


class TestB2:
    @pytest.mark.parametrize("radii", [[1.0, 1.0], [0.5, 1.5], [0.8, 1.2]])
    def test_returns_valid_tuple(self, radii):
        result = b2(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = b2([1.0, 0.5])
        a, b, c, *_ = result
        assert a == b == c


class TestZincblende:
    @pytest.mark.parametrize("radii", [[1.0, 0.5], [0.8, 0.8], [1.5, 0.3]])
    def test_returns_valid_tuple(self, radii):
        result = zincblende(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = zincblende([1.0, 0.5])
        a, b, c, *_ = result
        assert a == b == c


class TestB10:
    @pytest.mark.parametrize("radii", [[1.0, 0.5], [0.8, 1.2], [1.5, 0.3]])
    def test_returns_valid_tuple(self, radii):
        result = b10(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_tetragonal_symmetry(self):
        result = b10([1.0, 0.5])
        a, b, c, *_ = result
        assert a == b
        # c is different from a in litharge
        assert c != a


class TestStuffedWurtzite:
    @pytest.mark.parametrize("radii", [[0.7, 0.5, 1.0], [1.0, 0.8, 1.2], [0.5, 0.3, 0.6]])
    def test_returns_valid_tuple(self, radii):
        result = stuffed_wurtzite(radii)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, HEX_ANGLES)

    def test_hexagonal_symmetry(self):
        result = stuffed_wurtzite([0.7, 0.5, 1.0])
        a, b, _c, *_ = result
        assert a == b


# ---------------------------------------------------------------------------
# Tests for single-radius functions
# ---------------------------------------------------------------------------


class TestFcc:
    @pytest.mark.parametrize("radius", [0.5, 1.0, 1.5, 2.0])
    def test_returns_valid_tuple(self, radius):
        result = fcc(radius)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = fcc(1.0)
        a, b, c, *_ = result
        assert a == b == c


class TestBcc:
    @pytest.mark.parametrize("radius", [0.5, 1.0, 1.5, 2.0])
    def test_returns_valid_tuple(self, radius):
        result = bcc(radius)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = bcc(1.0)
        a, b, c, *_ = result
        assert a == b == c


class TestHcp:
    @pytest.mark.parametrize("radius", [0.5, 1.0, 1.5, 2.0])
    def test_returns_valid_tuple(self, radius):
        result = hcp(radius)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, HEX_ANGLES)

    def test_hexagonal_symmetry(self):
        result = hcp(1.0)
        a, b, _c, *_ = result
        assert a == b


class TestDiamond:
    @pytest.mark.parametrize("radius", [0.5, 1.0, 1.5, 2.0])
    def test_returns_valid_tuple(self, radius):
        result = diamond(radius)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_cubic_symmetry(self):
        result = diamond(1.0)
        a, b, c, *_ = result
        assert a == b == c


class TestBct:
    @pytest.mark.parametrize("radius", [0.5, 1.0, 1.5, 2.0])
    def test_returns_valid_tuple(self, radius):
        result = bct(radius)
        _assert_six_tuple(result)
        _assert_positive_lengths(result)
        _assert_angles(result, CUBIC_ANGLES)

    def test_tetragonal_symmetry(self):
        result = bct(1.0)
        a, b, c, *_ = result
        assert a == b
        assert c != a
