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
"""
from absl import app
from absl import flags
import enum
from fontTools.misc.bezierTools import splitCubic, splitCubicAtT
from functools import partial
from lxml import etree

from pathlib import Path

from picosvg.svg import SVG
from picosvg.geometric_types import Rect

def _line_pos(t, start, end):
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


def cubic_deriv_pos(t, p0, p1, p2, p3):
    return quad_pos(t, *cubic_deriv(p0, p1, p2, p3))


class FlagWarp:
  def __init__(self, minx, maxx):
    # cubic start, control, control, end
    # y-flipped for svg use
    self.warp = (
        (int(minx), 0),
        (int(minx + maxx / 3), -21),
        (int(minx + 2 * maxx / 3), 13),
        (int(maxx), -8)
    )
    self.minx = int(minx)
    self.maxx = int(maxx)
    self.miny = min(y for _, y in self.warp)
    self.maxy = max(y for _, y in self.warp)
    print(f"x [{self.minx}, {self.maxx}]")
    print(f"y [{self.miny}, {self.maxy}]")
    print("warp", self.warp)

  def _seg_ending_at(self, x):
    # Draw a line at pt.x through the warp; the pt on the warp is be the desired offset
    segments = splitCubic(*self.warp, x, False)

    # no kinks/overlaps; there should be <=2 curves when split
    # we are really just interested in there being only one unique start/end at pt.x
    assert len(segments) <= 2, f"Too many segments: {segments}"

    return segments[0]

  def _should_warp(self, pt):
    return pt[0] >= self.minx and pt[0] <= self.maxx

  def vec(self, pt):
    if not self._should_warp(pt):
        raise ValueError(f"{pt} outside [{self.minx}, {self.maxx}]")
    cubic_to_x = self._seg_ending_at(pt[0])
    # no x-movement, we just want the y-part of the last point
    result = (0, cubic_to_x[-1][1])
    assert result[1] >= self.miny and result[1] <= self.maxy
    return result

  def deriv(self, pt):
    if not self._should_warp(pt):
        return None
    return _cubic_deriv_pos(1, *self._seg_ending_at(pt[0]))


def _reduce_text(text):
    text = text.strip() if text else None
    return text if text else None


def _cubic_callback(num_splits, subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args):
  # Convert to cubic if needed
  if cmd.upper() == 'Z':
    # a line back to subpath start
    cmd = 'L'
    args = subpath_start

  if cmd.upper() == 'M':
    return ((cmd, args),)
  elif cmd == 'L':
    # line from curr_xy to args
    assert len(args) == 2
    end_xy = args

    # cubic ctl points 1/3rd and 2/3rds along
    cmd = "C"
    args = (
      *_line_pos(0.33, curr_xy, end_xy),
      *_line_pos(0.66, curr_xy, end_xy),
      *end_xy,
    )

  if cmd != 'C':
    raise ValueError(f"How do you cubic {cmd}, {args}")

  assert len(args) == 6

  split_t = tuple(t / num_splits for t in range(1, num_splits))
  cubic = (curr_xy, args[0:2], args[2:4], args[4:6])
  new_cmds = []
  for split_cubic in splitCubicAtT(*cubic, *split_t):
    assert len(split_cubic) == 4, split_cubic
    new_cmds.append(('C', split_cubic[1] + split_cubic[2] + split_cubic[3]))
  return tuple(new_cmds)


def _warp_callback(warp, subpath_start, curr_xy, cmd, args, prev_xy, prev_cmd, prev_args):
  args = list(args)
  new_args = []
  while(args):
    x, y = args.pop(0), args.pop(0)
    dx, dy = warp.vec((x, y))
    new_args.append(x + dx)
    new_args.append(y + dy)
    #print(f"({x:.2f}, {y:.2f}) => ({x+dx:.2f}, {y+dy:.2f}) delta ({dx:.2f}, {dy:.2f})")

  return ((cmd, tuple(new_args)),)


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
      max_corner = tuple(max(v1, v2) for v1, v2 in zip(max_corner, (box.x + box.w, box.y + box.h)))

    return Rect(*min_corner, *(maxv - minv for minv, maxv in zip(min_corner, max_corner)))


def _coordstr(c):
  return f"{c:.2f}"

def _dot(parent, at, radius=1):
  dot = etree.SubElement(parent, "circle")
  dot.attrib["fill"] = "cyan"
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

    warp = FlagWarp(box.x, box.x + box.w)
    warp_callback = partial(_warp_callback, warp)

    cubic_callback = partial(_cubic_callback, 4)

    
    for shape in svg.shapes():
      shape.explicit_lines(inplace=True)
      shape.walk(cubic_callback)
      shape.walk(warp_callback)
      shape.round_floats(2, inplace=True)

    tree = svg.toetree()

    # debug constructs
    for pt in warp.warp:
      _dot(tree, pt)

    for x in range(warp.minx, warp.maxx + 1, 4):
      p0 = (x, 0)
      p1 = (p0[0], warp.vec(p0)[1])
      _line(tree, p0, p1)

      p0 = (p0[0], box.y + box.h)
      p1 = (p1[0], p1[1] + box.y + box.h)
      _line(tree, p0, p1)

    # lxml really likes to retain whitespace
    for e in tree.iter("*"):
        e.text = _reduce_text(e.text)
        e.tail = _reduce_text(e.tail)

    print(etree.tostring(tree, pretty_print=True).decode("utf-8"))


if __name__ == "__main__":
    app.run(main)
