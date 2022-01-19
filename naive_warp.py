"""
Naive warp:
  pico the svg
  convert everything to cubic
  split the cubics (we'll just assume too coarse initially)
  warp start/end and intermediate points directly
    no fancy tricks with derivatives... my notes said we should have some

Note the lack of:
  curve simplification
  intelligent splitting (only when error introduced)

Usage:
  python naive_warp.py 1f1e6_1f1ec.svg
  python naive_warp.py noto-emoji/third_party/region-flags/svg/GB-WLS.svg --out_file GB-WLS-waved.svg
"""
import re
import sys
from absl import app
from absl import flags
from cu2qu import curve_to_quadratic as cubic_to_quad
import enum
from fontTools.misc.bezierTools import (
    approximateCubicArcLength,
    approximateQuadraticArcLength,
    splitCubic,
    splitCubicAtT,
    splitQuadraticAtT,
)
from functools import partial
from helpers import *
from gg_fit_curve import fit_cubics
from lxml import etree
from math import ceil, floor, sqrt, degrees, asin
from pathlib import Path
import pathops

from picosvg.geometric_types import Vector, Point, Rect
from picosvg.svg_meta import num_args, ntos, strip_ns
from picosvg.svg_transform import Affine2D
from picosvg.svg import _GRADIENT_CLASSES, SVGLinearGradient, SVGRadialGradient


DEFAULT_PRECISION = 1000
DEFAULT_FLATNESS = 1.0001
# The noto-emoji flags width/height aspect ratio is ~= 1.46. I found this empirically
# by measuring 126 / 86 pixels, centered horizontally and vertically within a 128 x 128
# PNG viewport. The values below are adjusted for 1000 x 1000 viewbox size.
DEFAULT_VIEWBOX_SIZE = 1000
DEFAULT_WIDTH = 984
DEFAULT_HEIGHT = 672
DEFAULT_RIGHT_MARGIN = 8
DEFAULT_TOP_MARGIN = 164
DEFAULT_BORDER_WIDTH = 1 / 32
# The aspect ratio of (the majority of) the input SVG flags
STD_ASPECT = 5 / 3

FLAGS = flags.FLAGS

flags.DEFINE_string("out_file", "-", "Output, - means stdout")
flags.DEFINE_bool("debug_info", False, "Add potentially useful debug marks")
flags.DEFINE_enum(
    "mode",
    "schneider_cubic",
    ["schneider_cubic", "quad", "dumb_cubic"],
    "How to generate final curves.",
)
flags.DEFINE_float(
    "precision", DEFAULT_PRECISION, "Default: 1/1000th of the viewbox diagonal"
)
flags.DEFINE_float(
    "flatness",
    DEFAULT_FLATNESS,
    "How closely a curve approximates a line: 1.0 exactly flat, > 1.0 flat enough",
)
flags.DEFINE_float("viewbox_size", DEFAULT_VIEWBOX_SIZE, "Viewbox width and height")
flags.DEFINE_float("width", DEFAULT_WIDTH, "Flag width")
flags.DEFINE_float("height", DEFAULT_HEIGHT, "Flag height")
flags.DEFINE_float("right_margin", DEFAULT_RIGHT_MARGIN, "Flag right margin")
flags.DEFINE_float("top_margin", DEFAULT_TOP_MARGIN, "Flag top margin")
flags.DEFINE_float(
    "border_size", DEFAULT_BORDER_WIDTH, "Relative border size as proportion of viewbox"
)


def line_pos(t, start, end):
    sx, sy = start
    ex, ey = end
    return (sx + t * (ex - sx), sy + t * (ey - sy))


# naive impl, essentially encoding something like:
# https://en.wikipedia.org/wiki/B%C3%A9zier_curve#/media/File:B%C3%A9zier_3_big.svg
def cubic_pos(t, p0, p1, p2, p3):

    q0 = line_pos(t, p0, p1)
    q1 = line_pos(t, p1, p2)
    q2 = line_pos(t, p2, p3)

    r0 = line_pos(t, q0, q1)
    r1 = line_pos(t, q1, q2)
    b = line_pos(t, r0, r1)

    return b


# https://pomax.github.io/bezierinfo/#derivatives
def cubic_deriv(p0, p1, p2, p3):
    # derivative is a quad
    return (
        Point(3 * (p1[0] - p0[0]), 3 * (p1[1] - p0[1])),
        Point(3 * (p2[0] - p1[0]), 3 * (p2[1] - p1[1])),
        Point(3 * (p3[0] - p2[0]), 3 * (p3[1] - p2[1])),
    )


def quad_pos(t, p0, p1, p2):
    q0 = line_pos(t, p0, p1)
    q1 = line_pos(t, p1, p2)
    b = line_pos(t, q0, q1)
    return b


def bez_pos(t, *points):
    if len(points) == 2:
        return line_pos(t, *points)
    elif len(points) == 3:
        return quad_pos(t, *points)
    elif len(points) == 4:
        return cubic_pos(t, *points)
    else:
        raise ValueError("Not enough points")


def cubic_deriv_pos(t, p0, p1, p2, p3):
    return quad_pos(t, *cubic_deriv(p0, p1, p2, p3))


def cubic_tangent(t, p0, p1, p2, p3) -> Vector:
    # Returns the unit vector defining the cubic bezier curve's direction at t
    if t == 0.0:
        tangent = Point(*p1) - Point(*p0)
    elif t == 1.0:
        tangent = Point(*p3) - Point(*p2)
    else:
        tangent = Vector(*cubic_deriv_pos(t, p0, p1, p2, p3))
    tangent = tangent.unit()
    if tangent is not None:
        return tangent
    return Vector()


def clamp(minv, maxv, v):
    if v < minv:
        return minv
    if v > maxv:
        return maxv
    return v


def quad_from_three_points(p0, p1, p2):
    # Construct a quadratic curve from three on-curve points: p0 and p2 being
    # the end-points, and p1 the on-curve point at t=0.5.
    # References:
    # https://pomax.github.io/bezierinfo/#abc
    # https://pomax.github.io/bezierinfo/#pointcurves

    B = p1

    u = 0.5  # projection ratio for t=0.5 and order n=2

    # given the above ratio, C is always half-way between p0 and p2
    C = (u * p0[0] + (1 - u) * p2[0], u * p0[1] + (1 - u) * p2[1])

    s = 1.0  # ABC ratio for t=0.5 and order n=2

    # compute the quadratic off-curve point A given B, C and ABC ratio
    A = (B[0] + (B[0] - C[0]) / s, B[1] + (B[1] - C[1]) / s)

    return (p0, A, p2)


class FlagWarp:
    def __init__(self, box, viewbox_size):
        # cubic start, control, control, end
        # based on Noto Emoji waveflag.c warp
        y_scale = viewbox_size / 128
        self.minx = floor(box.x)
        self.maxx = ceil(box.x + box.w)
        self.warp = (
            (self.minx, 4 * y_scale),
            (self.minx + box.w / 3, -17 * y_scale),
            (self.minx + 2 * box.w / 3, 17 * y_scale),
            (self.maxx, -4 * y_scale),
        )
        self.miny = min(y for _, y in self.warp)
        self.maxy = max(y for _, y in self.warp)

    def _seg_ending_at(self, x):
        x = clamp(self.minx, self.maxx, x)

        # Draw a line at pt.x through the warp; the pt on the warp is be the desired offset
        segments = splitCubic(*self.warp, x, False)

        # no kinks/overlaps; there should be <=2 curves when split
        # we are really just interested in there being only one unique start/end at pt.x
        assert len(segments) <= 2, f"Too many segments: {segments}"

        return segments[0]

    def vec(self, pt):
        cubic_to_x = self._seg_ending_at(pt[0])
        # no x-movement, we just want the y-part of the last point
        result = (0, cubic_to_x[-1][1])
        assert result[1] >= self.miny and result[1] <= self.maxy
        return result

    def skewY(self, x: float) -> Affine2D:
        # Get linear transform that approximates the flag's (non-linear) warp around
        # the given x coordinate ("Jacobian matrix"?).
        if x <= self.minx:
            tangent = cubic_tangent(0.0, *self.warp)
        else:
            tangent = cubic_tangent(1.0, *self._seg_ending_at(x))
        angle = degrees(asin(tangent.y))
        return Affine2D.fromstring(f"translate({x}) skewY({angle}) translate({-x})")


def close_open_subpaths(path):
    def callback(subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args):
        if cmd.upper() == "M" and prev_cmd is not None and prev_cmd.upper() != "Z":
            return (("Z", ()), (cmd, args))
        return ((cmd, args),)

    path.walk(callback)
    if not path.d.endswith(("Z", "z")):
        path.end()


def _cubic_callback(subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args):
    if cmd.upper() == "M":
        return ((cmd, args),)

    # Convert to cubic if needed
    if cmd.upper() == "Z":
        # a line back to subpath start ... unless we are there already
        if curr_xy == subpath_start:
            return ()
        cmd = "L"
        args = subpath_start

    if cmd == "L":
        # line from curr_xy to args
        assert len(args) == 2
        end_xy = args

        # cubic ctl points 1/3rd and 2/3rds along
        cmd = "C"
        args = (
            *line_pos(0.33, curr_xy, end_xy),
            *line_pos(0.66, curr_xy, end_xy),
            *end_xy,
        )

    if cmd == "Q":
        assert len(args) == 4
        cmd = "C"
        p0 = Point(*curr_xy)
        p1 = Point(*args[0:2])
        p2 = Point(*args[2:4])
        args = (*(p0 + 2 / 3 * (p1 - p0)), *(p2 + 2 / 3 * (p1 - p2)), *p2)

    if cmd != "C":
        raise ValueError(f"How do you cubic {cmd}, {args}")

    assert len(args) == 6

    return (("C", args),)


def _quadratic_callback(
    max_err, subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args
):
    if cmd.upper() == "M":
        return ((cmd, args),)

    # Convert to quadratic if needed
    if cmd.upper() == "Z":
        # a line back to subpath start ... unless we are there already
        if curr_xy == subpath_start:
            return ()
        cmd = "L"
        args = subpath_start

    if cmd == "L":
        # line from curr_xy to args
        assert len(args) == 2
        end_xy = args

        # quad ctl point 1/2 along
        cmd = "Q"
        args = (
            *line_pos(0.5, curr_xy, end_xy),
            *end_xy,
        )

    if cmd == "C":
        cubic_points = [curr_xy] + [
            tuple(args[i : i + 2]) for i in range(0, len(args), 2)
        ]
        quad_points = cubic_to_quad(cubic_points, max_err)
        quad_segments = []
        for (control_pt, end_pt) in pathops.decompose_quadratic_segment(
            tuple(quad_points[1:])
        ):
            quad_segments.append(("Q", (*control_pt, *end_pt)))
        return tuple(quad_segments)

    if cmd != "Q":
        raise ValueError(f"How do you quad {cmd}, {args}")

    assert len(args) == 4

    return (("Q", args),)


def _dist(p0, p1):
    return sqrt(abs(p0[0] - p1[0]) ** 2 + abs(p0[1] - p1[1]))


def points(curve, num_seg=100):
    # num_seq+1 coords on curve, equidistent in t (1/(num_seg) apart)
    # basically http://pomax.github.io/bezierjs/#getLUT
    return tuple(bez_pos(step / num_seg, *curve) for step in range(0, num_seg + 1))


def _max_pt_err(view_box, precision=DEFAULT_PRECISION):
    # TODO is this remotely reasonable?
    return _dist((0, 0), (view_box.w, view_box.h)) / precision


def _quad_max_err(view_box):
    min_length = min(view_box.w, view_box.h)
    return min_length / 1000


def _ref_pts(view_box, curve):
    # we award you one segment per 1000th of viewbox diagonal
    # 100th didn't seem to be enough
    if len(curve) == 4:
        curve_len = approximateCubicArcLength(*curve)
    elif len(curve) == 3:
        curve_len = approximateQuadraticArcLength(*curve)
    else:
        raise AssertionError(len(curve))
    diagonal = _dist((0, 0), (view_box.w, view_box.h))
    num_seg = max(ceil(1000 * curve_len / diagonal), 10)
    return points(curve, num_seg)


def _nearest(p0, pts):
    result = (p0, pts[0])
    for p1 in pts[1:]:
        if _dist(p0, p1) < _dist(*result):
            result = (p0, p1)
    return result


def _bez_nearest_range_t(curve, curr_t=0, num_seg=100, min_t=0, max_t=1):
    min_t = clamp(0, 1, min_t)
    max_t = clamp(0, 1, max_t)

    curr_t = clamp(min_t, max_t, curr_t)
    curr_dist = _dist(bez_pos(curr_t, *curve))
    for step in range(0, num_seg + 1):
        try_t = clamp(min_t, max_t, min_t + step * (max_t - min_t) / num_seg)
        try_dist = _dist(bez_pos(curr_t, *curve), bez_pos(try_t, *curve))
        if try_dist < curr_dist:
            curr_t = try_t
            curr_dist = try_dist
    return curr_t


def _bez_nearest_t(curve, pt):
    # TODO: kurbo has a real implementation but this will do for experimentation
    # inspired by http://pomax.github.io/bezierjs/#project

    # try over full range then try to refine answer
    curr_t = _bez_nearest_range_t(curve, curr_t=0, num_seg=100)

    # we could repeat this as long as answer improves sufficiently?
    curr_t = _bez_nearest_range_t(
        curve, curr_t=0, num_seg=100, min_t=curr_t - 1 / 100, max_t=curr_t + 1 / 100
    )

    return curr_t


def _bez_nearest(curve, pt):
    return bez_pos(_bez_nearest_t(curve, pt), *curve)


def _warp_pt(warp, pt):
    dx, dy = warp.vec(pt)
    return (pt[0] + dx, pt[1] + dy)


def _unwarp_pt(warp, pt):
    dx, dy = warp.vec(pt)
    return (pt[0] - dx, pt[1] - dy)


def _warp_cubic(warp, cubic):
    assert len(cubic) == 4
    return tuple(_warp_pt(warp, pt) for pt in cubic)


def _warp_quad(warp, quad):
    assert len(quad) == 3
    quad_pts = (quad[0], quad_pos(0.5, *quad), quad[2])
    return quad_from_three_points(*(_warp_pt(warp, p) for p in quad_pts))


class PathWarp:
    def __init__(
        self, view_box, warp, cmd, warp_curve_fn, split_at_t_fn, precision, flatness
    ):
        self._view_box = view_box
        self._warp = warp
        self._cmd = cmd
        self._warp_curve_fn = warp_curve_fn
        self._split_at_t_fn = split_at_t_fn
        self._max_pt_err = _max_pt_err(view_box, precision)
        self._flatness = flatness

    def warp_callback(
        self, subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args
    ):
        if cmd.upper() == "M":
            return ((cmd, _warp_pt(self._warp, args)),)

        assert cmd == self._cmd, f"{cmd} != {self._cmd}"
        assert len(args) == num_args(cmd)

        # build reference points (# based on arc length)
        # unwarp: current was already warped, but we'll warp again when processing
        initial_curve = (_unwarp_pt(self._warp, curr_xy),) + tuple(
            args[i * 2 : (i + 1) * 2] for i in range(int(len(args) / 2))
        )
        ref_pts = tuple(
            _warp_pt(self._warp, pt) for pt in _ref_pts(self._view_box, initial_curve)
        )

        frontier = [initial_curve]
        acceptable = []
        while frontier:
            # take one down...
            curve = frontier.pop(0)
            # warp it around...
            warped_curve = self._warp_curve_fn(self._warp, curve)

            # ...good enough?
            warp_refs = _ref_pts(self._view_box, warped_curve)
            ok = all(
                _dist(*_nearest(pt, ref_pts)) <= self._max_pt_err for pt in warp_refs
            )

            if ok:
                acceptable.append(warped_curve)
            else:
                # Cut the pre-warp curve in half and try again
                # Intuition: there are better guesses than half
                splits = list(self._split_at_t_fn(*curve, 0.5))
                assert len(splits) > 1
                frontier = splits + frontier
            # print(
            #     f"frontier {len(frontier)} accepted {len(acceptable)}; {len(warp_refs)} for warp"
            # )

        return tuple(
            ("L", q[-1])
            if _is_almost_line(q, self._flatness)
            else (cmd, sum(q[1:], ()))
            for q in acceptable
        )


def _quad_path_warp(view_box, warp, precision, flatness):
    return PathWarp(
        view_box, warp, "Q", _warp_quad, splitQuadraticAtT, precision, flatness
    )


def _cubic_path_warp(view_box, warp, precision, flatness):
    return PathWarp(
        view_box, warp, "C", _warp_cubic, splitCubicAtT, precision, flatness
    )


def _is_almost_line(curve, flatness=1.0):
    # Returns True if the bezier curve is equivalent to a line.
    # A 'flat' bezier curve is one such that the sum of the distances between
    # consecutive control points equals the distance from start to end points.
    # That's because if a control point isn't on the line from start to end then
    # the length would exceed the direct path start => end.
    # The 'flatness' factor of 1.0 means exactly flat, anything greater than 1.0
    # proportionately means flat "enough". Less than 1.0 means never flat (i.e.
    # keep all flat curves as curves, FWIW).
    points = [Point(*p) for p in curve]
    length = 0
    for i in range(len(curve) - 1):
        length += (points[i + 1] - points[i]).norm()
    max_length = (points[-1] - points[0]).norm()
    return length <= flatness * max_length


class FitCubicPathWarp:
    def __init__(self, view_box, warp, precision, flatness):
        self._view_box = view_box
        self._warp = warp
        # fitCurves.py uses squared distances
        self._max_err = _max_pt_err(view_box, precision) ** 2
        self._flatness = flatness

    def warp_callback(
        self, subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args
    ):
        if cmd.upper() == "M":
            return ((cmd, _warp_pt(self._warp, args)),)

        assert cmd == "C", f"{cmd} != {self._cmd}"
        assert len(args) == num_args(cmd)

        # build reference points (# based on arc length)
        # unwarp: current was already warped, but we'll warp again when processing
        initial_curve = (_unwarp_pt(self._warp, curr_xy),) + tuple(
            args[i * 2 : (i + 1) * 2] for i in range(int(len(args) / 2))
        )
        warped_pts = tuple(
            _warp_pt(self._warp, pt) for pt in _ref_pts(self._view_box, initial_curve)
        )

        # sampled points are evenly spaced over time, we want them to keep the same
        # parametrization after warping. E.g. the point at t=0.5 before warping will
        # still be located at t=0.5 on the warped curve.
        warped_curves = fit_cubics(warped_pts, self._max_err, uniform_parameters=True)

        return tuple(
            ("L", c[-1])
            if _is_almost_line(c, self._flatness)
            else (cmd, sum((tuple(p) for p in c[1:]), ()))
            for c in warped_curves
        )


def _bbox(boxes):
    min_corner = (9999, 9999)
    max_corner = (-9999, -9999)
    for box in boxes:
        min_corner = tuple(min(v1, v2) for v1, v2 in zip(min_corner, (box.x, box.y)))
        max_corner = tuple(
            max(v1, v2) for v1, v2 in zip(max_corner, (box.x + box.w, box.y + box.h))
        )

    return Rect(
        *min_corner, *(maxv - minv for minv, maxv in zip(min_corner, max_corner))
    )


def _coordstr(c):
    return f"{c:.2f}"


def _dot(parent, at, radius=1, color="cyan"):
    dot = etree.SubElement(parent, "circle")
    dot.attrib["fill"] = color
    dot.attrib["cx"] = _coordstr(at[0])
    dot.attrib["cy"] = _coordstr(at[1])
    dot.attrib["r"] = _coordstr(radius)


def _line(parent, p0, p1):
    line = etree.SubElement(parent, "line")
    line.attrib["stroke"] = "orange"
    line.attrib["x1"] = _coordstr(p0[0])
    line.attrib["y1"] = _coordstr(p0[1])
    line.attrib["x2"] = _coordstr(p1[0])
    line.attrib["y2"] = _coordstr(p1[1])


def _mid_point(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return Point(x1, y1) + (Point(x2, y2) - Point(x1, y1)) * 0.5


def _gradient_mid_point(gradient):
    if isinstance(gradient, SVGLinearGradient):
        return _mid_point((gradient.x1, gradient.y1), (gradient.x2, gradient.y2))
    elif isinstance(gradient, SVGRadialGradient):
        return _mid_point((gradient.cx, gradient.cy), (gradient.fx, gradient.fy))
    else:
        raise TypeError(type(gradient))


def _transform_gradients(self, warp):
    for el in self._select_gradients():
        gradient = _GRADIENT_CLASSES[strip_ns(el.tag)].from_element(el, self.view_box())

        mid_pt = _gradient_mid_point(gradient)
        mid_pt = gradient.gradientTransform.map_point(mid_pt)
        translate = Affine2D.identity().translate(*warp.vec(mid_pt))
        mid_pt = translate.map_point(mid_pt)

        gradient_transform = Affine2D.compose_ltr(
            (gradient.gradientTransform, translate, warp.skewY(mid_pt.x))
        )

        if gradient_transform != Affine2D.identity():
            el.attrib["gradientTransform"] = gradient_transform.tostring()
        elif "gradientTransform" in el.attrib:
            del el.attrib["gradientTransform"]
    return self


# TODO: move to picosvg proper?
# First need to fix https://github.com/googlefonts/picosvg/issues/110
def _picosvg_transform(self, affine):
    for idx, (el, (shape,)) in enumerate(self._elements()):
        self.elements[idx] = (el, (shape.apply_transform(affine),))

    for el in self._select_gradients():
        gradient = _GRADIENT_CLASSES[strip_ns(el.tag)].from_element(el, self.view_box())
        gradient_transform = Affine2D.compose_ltr((gradient.gradientTransform, affine))

        if gradient_transform != Affine2D.identity():
            el.attrib["gradientTransform"] = gradient_transform.tostring()
        elif "gradientTransform" in el.attrib:
            del el.attrib["gradientTransform"]
    return self


SVG_UNITS_RE = re.compile("(?:em|ex|px|pt|pc|cm|mm|in|%)$")


# TODO: move to picosvg?
def apply_viewbox_preserve_aspect_ratio(svg):
    """If viewport != viewBox apply the resulting transform and remove viewBox.
    Takes 'preserveAspectRatio' into account.
    E.g. The Qatar flag (QA.svg) needs this treatment.
    """
    svg_root = svg.svg_root
    width = svg_root.attrib.get("width")
    height = svg_root.attrib.get("height")
    if width is not None and height is not None and "viewBox" in svg_root.attrib:
        # ignore absolute length units; we're only interested in the relative size
        # of viewport vs viewbox here
        width = SVG_UNITS_RE.sub("", width)
        height = SVG_UNITS_RE.sub("", height)
        viewport = Rect(0, 0, float(width), float(height))
        viewbox = svg.view_box()
        if viewport != viewbox:
            transform = Affine2D.rect_to_rect(
                viewbox,
                viewport,
                svg_root.attrib.get("preserveAspectRatio", "xMidYMid"),
            )
            _picosvg_transform(svg, transform)
            del svg_root.attrib["viewBox"]
            if "preserveAspectRatio" in svg_root.attrib:
                del svg_root.attrib["preserveAspectRatio"]


def _x_aspect(v, aspect, size):
    # scale abscissa around mid-point if aspect < 1.0 (i.e. width narrower than height)
    mid = size / 2
    return v if aspect >= 1.0 else (v - mid) * aspect + mid


def _y_aspect(v, aspect, size):
    # scale ordinate around mid-point if aspect > 1.0 (i.e. height taller than width)
    mid = size / 2
    return v if aspect <= 1.0 else (v - mid) / aspect + mid


def normalize_flag_aspect(svg, viewbox_size, width, height, right_margin, top_margin):
    apply_viewbox_preserve_aspect_ratio(svg)

    current_box = _bbox(tuple(s.bounding_box() for s in svg.shapes()))

    # Try to keep overall proportions for the flags that are considerably
    # narrower or wider than the standard aspect ratio.
    aspect = current_box.w / current_box.h
    aspect /= STD_ASPECT
    aspect = sqrt(aspect)  # Discount the effect
    if 0.9 <= aspect <= 1.1:
        aspect = 1.0
    else:
        print("Non-standard aspect ratio:", aspect)

    xmin = _x_aspect(right_margin, aspect, viewbox_size)
    ymin = _y_aspect(top_margin, aspect, viewbox_size)
    xmax = _x_aspect(right_margin + width, aspect, viewbox_size)
    ymax = _y_aspect(top_margin + height, aspect, viewbox_size)
    new_box = Rect(xmin, ymin, xmax - xmin, ymax - ymin)

    affine = Affine2D.rect_to_rect(current_box, new_box)

    _picosvg_transform(svg, affine)

    square_viewbox = Rect(0, 0, viewbox_size, viewbox_size)
    svg.svg_root.attrib["viewBox"] = " ".join(ntos(v) for v in square_viewbox)
    for attr_name in ("width", "height"):
        if attr_name in svg.svg_root.attrib:
            del svg.svg_root.attrib[attr_name]


def make_border_svg(viewbox_size, border_size, box):
    stroke_width = border_size * viewbox_size * 2
    # this is to make horizontal strokes a bit thicker to compensate for the non-linear
    # warp that makes them look thinner than the vertical ones
    padding = stroke_width * 0.08
    svg = SVG.fromstring(
        f"""<svg version="1.1" xmlns="http://www.w3.org/2000/svg"
                     viewBox="0 0 {viewbox_size} {viewbox_size}">
              <defs>
                <clipPath id="clip">
                  <rect x="{box.x}" y="{box.y}"
                        width="{box.w}" height="{box.h}"/>
                </clipPath>
              </defs>
              <rect fill="none" stroke="#1A1A1A" stroke-width="{stroke_width}"
                    stroke-opacity="0.2" x="{box.x}" y="{box.y + padding / 2}"
                    width="{box.w}" height="{box.h - padding}" clip-path="url(#clip)"/>
            </svg>"""
    ).topicosvg(inplace=True)
    return svg


def main(argv):
    svg = load_svg(argv).topicosvg().clip_to_viewbox(inplace=True)

    normalize_flag_aspect(
        svg,
        FLAGS.viewbox_size,
        FLAGS.width,
        FLAGS.height,
        FLAGS.right_margin,
        FLAGS.top_margin,
    )

    box = _bbox(tuple(s.bounding_box() for s in svg.shapes()))
    warp = FlagWarp(box, viewbox_size=FLAGS.viewbox_size)
    precision = FLAGS.precision
    flatness = FLAGS.flatness

    print(f"{FLAGS.mode} mode...", file=sys.stderr)

    if FLAGS.mode == "schneider_cubic":
        prep_callback = _cubic_callback
        path_warp = FitCubicPathWarp(svg.view_box(), warp, precision, flatness)
    elif FLAGS.mode == "quad":
        prep_callback = partial(_quadratic_callback, _quad_max_err(svg.view_box()))
        path_warp = _quad_path_warp(svg.view_box(), warp, precision, flatness)
    elif FLAGS.mode == "dumb_cubic":
        prep_callback = _cubic_callback
        path_warp = _cubic_path_warp(svg.view_box(), warp, precision, flatness)
    else:
        raise ValueError("Specify a valid --mode")

    if FLAGS.border_size:
        border = make_border_svg(FLAGS.viewbox_size, FLAGS.border_size, box)
        svg.append_to("/svg:svg", border.xpath_one("//svg:path"))

    for shape in svg.shapes():
        shape.explicit_lines(inplace=True)
        shape.arcs_to_cubics(inplace=True)
        shape.expand_shorthand(inplace=True)
        close_open_subpaths(shape)
        shape.walk(prep_callback)
        shape.walk(path_warp.warp_callback)
        shape.round_floats(2, inplace=True)

    _transform_gradients(svg, warp)

    tree = svg.toetree()

    # debug constructs
    if FLAGS.debug_info:
        for pt in warp.warp:
            _dot(tree, pt)

        for x in range(warp.minx, warp.maxx + 1, 4):
            p0 = (x, 0)
            p1 = (p0[0], warp.vec(p0)[1])
            _line(tree, p0, p1)

            p0 = (p0[0], box.y + box.h)
            p1 = (p1[0], p1[1] + box.y + box.h)
            _line(tree, p0, p1)

        _line(tree, (0, box.y), (138, box.y))
        _line(tree, (0, box.y + box.h), (138, box.y + box.h))

        for shape in svg.shapes():
            for cmd, args in shape.as_cmd_seq():
                for i in range(0, len(args), 2):
                    color = "gray"
                    if cmd.upper() == "M":
                        color = "darkgreen"
                    elif i == len(args) - 1:
                        color = "black"
                    pt = args[i : i + 2]
                    _dot(tree, pt, radius=1.5, color=color)

                    color = "blue"
                    wv = warp.vec(args[i : i + 2])
                    pt = (pt[0] - wv[0], pt[1] - wv[1])
                    # _dot(tree, pt, radius=1.5, color=color)
        _dot(tree, (0, 0), radius=2.5, color="purple")

    reduce_whitespace(tree)
    write_xml(FLAGS.out_file, tree)


if __name__ == "__main__":
    app.run(main)
