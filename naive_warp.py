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
import sys
from absl import app
from absl import flags
from cu2qu import curve_to_quadratic as cubic_to_quad
import enum
from fontTools.misc.bezierTools import (
    calcCubicArcLength,
    calcQuadraticArcLength,
    splitCubic,
    splitCubicAtT,
    splitQuadraticAtT,
)
from functools import partial
from gg_fit_curve import fit_cubics
from lxml import etree
from math import ceil, floor, sqrt
from pathlib import Path
import pathops

from picosvg.geometric_types import Rect
from picosvg.svg import SVG
from picosvg.svg_meta import num_args


DEFAULT_PRECISION = 200

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
    "precision", DEFAULT_PRECISION, "Default: 1/200th of the viewbox diagonal"
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
        Point(3 * (p1.x - p0.x), 3 * (p1.y - p0.y)),
        Point(3 * (p2.x - p1.x), 3 * (p2.y - p1.y)),
        Point(3 * (p3.x - p2.x), 3 * (p3.y - p2.y)),
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
    def __init__(self, box):
        # cubic start, control, control, end
        # based on Noto Emoji waveflag.c warp
        y_scale = box.w / 128
        self.minx = floor(box.x)
        self.maxx = ceil(box.x + box.w)
        self.warp = (
            (self.minx, 0),
            (self.minx + self.maxx / 3, -21 * y_scale),
            (self.minx + 2 * self.maxx / 3, 13 * y_scale),
            (self.maxx, -8 * y_scale),
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

    def deriv(self, pt):
        return _cubic_deriv_pos(1, *self._seg_ending_at(pt[0]))


def _reduce_text(text):
    text = text.strip() if text else None
    return text if text else None


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
        curve_len = calcCubicArcLength(*curve)
    elif len(curve) == 3:
        curve_len = calcQuadraticArcLength(*curve)
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
    def __init__(self, view_box, warp, cmd, warp_curve_fn, split_at_t_fn, precision):
        self._view_box = view_box
        self._warp = warp
        self._cmd = cmd
        self._warp_curve_fn = warp_curve_fn
        self._split_at_t_fn = split_at_t_fn
        self._max_pt_err = _max_pt_err(view_box, precision)

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

        return tuple((cmd, sum(q[1:], ())) for q in acceptable)


def _quad_path_warp(view_box, warp, precision):
    return PathWarp(view_box, warp, "Q", _warp_quad, splitQuadraticAtT, precision)


def _cubic_path_warp(view_box, warp, precision):
    return PathWarp(view_box, warp, "C", _warp_cubic, splitCubicAtT, precision)


class FitCubicPathWarp:
    def __init__(self, view_box, warp, precision):
        self._view_box = view_box
        self._warp = warp
        # fitCurves.py uses squared distances
        self._max_err = _max_pt_err(view_box, precision) ** 2

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

        return tuple((cmd, sum((tuple(p) for p in c[1:]), ())) for c in warped_curves)


def _load_svg(argv):
    assert len(argv) <= 2
    try:
        input_file = argv[1]
    except IndexError:
        input_file = None

    if input_file:
        svg = SVG.parse(input_file)
    else:
        svg = SVG.fromstring(sys.stdin.read())
    return svg.topicosvg()


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


def main(argv):
    svg = _load_svg(argv)

    box = _bbox(tuple(s.bounding_box() for s in svg.shapes()))

    warp = FlagWarp(box)
    precision = FLAGS.precision

    print(f"{FLAGS.mode} mode...")

    if FLAGS.mode == "schneider_cubic":
        prep_callback = _cubic_callback
        path_warp = FitCubicPathWarp(svg.view_box(), warp, precision)
    elif FLAGS.mode == "quad":
        prep_callback = partial(_quadratic_callback, _quad_max_err(svg.view_box()))
        path_warp = _quad_path_warp(svg.view_box(), warp, precision)
    elif FLAGS.mode == "dumb_cubic":
        prep_callback = _cubic_callback
        path_warp = _cubic_path_warp(svg.view_box(), warp, precision)
    else:
        raise ValueError("Specify a valid --mode")

    for shape in svg.shapes():
        shape.explicit_lines(inplace=True)
        shape.walk(prep_callback)
        shape.walk(path_warp.warp_callback)
        shape.round_floats(2, inplace=True)

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

    # lxml really likes to retain whitespace
    for e in tree.iter("*"):
        e.text = _reduce_text(e.text)
        e.tail = _reduce_text(e.tail)

    out_content = etree.tostring(tree, pretty_print=True).decode("utf-8")
    if FLAGS.out_file == "-":
        print(out_content)
    else:
        with open(FLAGS.out_file, "w") as f:
            f.write(out_content)


if __name__ == "__main__":
    app.run(main)
