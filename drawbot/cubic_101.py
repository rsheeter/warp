from typing import NamedTuple
from fontTools.misc.bezierTools import splitCubic, splitCubicAtT

size(275, 200)
fontSize(6)

class Point(NamedTuple):
    x: float
    y: float

    def magnitude(self):
        return sqrt(self.x * self.x + self.y * self.y)
    
    def add(self, other):
        return Point(self.x + other.x, self.y + other.y)

    def of_magnitude(self, new_mag):
        # 0? hopefully not today :D
        mag = self.magnitude()
        return Point(new_mag * self.x / mag, new_mag * self.y / mag)   


def dot(at, radius=2):
    oval(at.x-radius/2, at.y-radius/2, radius, radius)


def line_pos(t, start, end):
    return Point(
        start.x + t * (end.x - start.x),
        start.y + t * (end.y - start.y),
    )

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


def bezOfCubics(cubics):
    bez = BezierPath()
    bez.moveTo(cubics[0][0])
    for cubic in cubics:
        bez.curveTo(*cubic[1:])
    return bez

def points(cubic, num_seg):
    # num_seq+1 coords on curve, equidistent in t
    return tuple(
        cubic_pos(t/num_seg, *cubic)
        for t in range(0, num_seg + 1)
    )

def draw_points(points, radius=1):
    for pt in points:
        dot(pt, radius)
      
# courtesy of @anthrotype, lightly modified
def cubic_to_hermite(cubic):
    # u0 and u1 are the tagents at each endpoint
    # courtesy of first Bezier derivative
    u0, _, u1 = cubic_deriv(*cubic)

    # endpt, tan, endpoint, tan
    return (cubic[0], u0, cubic[-1], u1)


# courtesy of @anthrotype, lightly modified
def hermite_to_cubic(hermite):
    # p1 and p2 are the Bezier control points, located at 1/3 distance of the
    # tangent's magnitude from each endpoint.
    # See "express cubic Hermite interpolation in terms of cubic Bézier curves" in
    # https://en.wikipedia.org/wiki/Cubic_Hermite_spline
    def ctl_pt(endpt, tanpt, mul):
        return Point(*(e + mul * t / 3 for e, t in zip(endpt, tanpt)))

    return (
        hermite[0],
        ctl_pt(*hermite[0:2], 1),
        ctl_pt(*hermite[3:1:-1], -1),
        hermite[3],
    )

c_test = (Point(0, 0), Point(1,0), Point(2,0), Point(3,0))
print("c_test", c_test)
hermite = cubic_to_hermite(c_test)
print("hermite", hermite)
c_again = hermite_to_cubic(hermite)
print("c_again", c_again)
assert c_test == c_again

def visualizeWarp(initial_cubics, warp_cubic_fn, warp_pt_fn, num_seg=10, draw_ctls=False):
    warped_cubics = []
    for cubic in initial_cubics:
        split_cubics = splitCubicAtT(*cubic, *(t/num_seg for t in range(0, num_seg+1)))
        for split_cubic in split_cubics:
            split_cubic = tuple(Point(*t) for t in split_cubic)
            # naive; just move according to shift
            # in the end expecting control point movement to differ from start/end point
            warped_cubic = warp_cubic_fn(split_cubic)
            warped_cubics.append(warped_cubic)
    warped_cubics = tuple(warped_cubics)

    fill(None)
    stroke(0.5, 0.5, 0.5)
    strokeWidth(0.25)
    for cubic in initial_cubics:
        drawPath(bezOfCubics((cubic,)))

    stroke(0.75, 0.5, 0.75)
    strokeWidth(0.4)
    for warp_cubic in warped_cubics:
        drawPath(bezOfCubics((warp_cubic,)))
        if draw_ctls:
            line(*warp_cubic[0:2])
            line(*warp_cubic[2:4])

    for cubic in initial_cubics:
        pts = tuple(warp_pt_fn(p) for p in points(cubic, 100))
        for i in range(1, len(pts)):
            line(pts[i-1], pts[i])
        for pt in pts:
            dot(pt, radius=1)

def fancyWarpPt(pt):
    mag_scale = 10
    val_scale = .25
    x,y = pt
    new_pt = Point(
        round(x + mag_scale * sin(val_scale * y), 2),
        round(y + mag_scale * sin(val_scale * x), 2)
    )
    return new_pt

def fancyWarpViaHermite(cubic):
    hermite = cubic_to_hermite(cubic)
    warp_hermite = tuple(fancyWarpPt(pt) for pt in hermite)
    warp_cubic = hermite_to_cubic(warp_hermite)
    return warp_cubic

##################################################
# simple hermite warp test
##################################################
simple_cubics = (
    (Point(50,75), Point(100,75), Point(150,75), Point(200,75)),
)

visualizeWarp(
    simple_cubics,
    lambda cubic: tuple(fancyWarpPt(pt) for pt in cubic),
    fancyWarpPt,
    num_seg=1,
    draw_ctls=True,
)

##################################################
# fancy warp test
##################################################
newPage()

# start with a series of linear cubics
cubics = []
for v in range(25, 200, 25):
    # horizontal
    cubics.append((Point(25, v), Point(75, v), Point(125, v), Point(175, v)))

    # vertical
    cubics.append((Point(v, 25), Point(v, 75), Point(v, 125), Point(v, 175)))

visualizeWarp(
    cubics,
    lambda cubic: tuple(fancyWarpPt(pt) for pt in cubic),
    fancyWarpPt
)

##################################################
# fancy hermite warp test
##################################################
newPage()

visualizeWarp(cubics, fancyWarpViaHermite, fancyWarpPt)

##################################################
# Points equidistant over t test
##################################################
newPage()

cubic = (
    Point(100,195), 
    Point(10,100),
    Point(180,20),
    Point(150,25)
)

fill(None)
stroke(0.5, 0.5, 0.5)
strokeWidth(0.25)
drawPath(bezOfCubics((cubic,)))

stroke(None)
fill(0, 0, 0.75)
draw_points(points(cubic, 16), radius=1.5)

fill(None)
stroke(0.5, 0.5, 0.5)
strokeWidth(0.2)
draw_points(cubic, 2.5)

##################################################
# Naive warp test
##################################################
newPage()

# from waveflag.c
top = 21
bot = 128-top
B = 21
C = 4


waveflag_mesh = (
  Point(  1, top+C),
  Point( 43, top-B+C),
  Point( 85, top+B-C),
  Point(127, top-C),
  Point(127, bot-C),
  Point( 85, bot+B-C),
  Point( 43, bot-B+C),
  Point(  1, bot+C),
)

oncurve = {0, 3, 4, 7}
offcurve = {1, 2, 5, 6}

# draw the basic shape and control points
translate(x=2, y=70)
text("Basic Shape", (0, 116))
#	cairo_line_to(cr,   M(0));
#	cairo_curve_to(cr,  M(1), M(2), M(3));
#	cairo_line_to(cr,   M(4));
#	cairo_curve_to(cr,  M(5), M(6), M(7));
#	cairo_close_path (cr);
path = BezierPath()
path.moveTo(waveflag_mesh[0])
path.curveTo(*waveflag_mesh[1:4])
path.lineTo(waveflag_mesh[4])
path.curveTo(*waveflag_mesh[5:8])
path.closePath()

fill(None)
stroke(1, 0, 0)
drawPath(path)

stroke(0.8, 0.8, 0.8)
strokeWidth(0.4)
line(waveflag_mesh[0], waveflag_mesh[1])
line(waveflag_mesh[1], waveflag_mesh[2])
line(waveflag_mesh[2], waveflag_mesh[3])

line(waveflag_mesh[4], waveflag_mesh[5])
line(waveflag_mesh[5], waveflag_mesh[6])
line(waveflag_mesh[6], waveflag_mesh[7])

# Oncurve points
fill(0, 0, 0)
stroke(None)
for idx in oncurve:
    dot(waveflag_mesh[idx])

# Offcurve points
fill(0.5, 0.5, 0.5)
stroke(None)
for idx in offcurve:
    pt = waveflag_mesh[idx]
    oval(pt.x-1, pt.y-1, 2, 2)    

# points along curve

mesh_deriv = cubic_deriv(*waveflag_mesh[0:4])
line(mesh_deriv[0], mesh_deriv[1])
line(mesh_deriv[1], mesh_deriv[2])

stroke(None)
num_pts = 20
for t in range(0, num_pts + 1):
    pt = cubic_pos(t/num_pts, *waveflag_mesh[0:4])
    fill(0, 0, 0.75)
    dot(pt, radius=1.5)
    
    pt = quad_pos(t, *mesh_deriv)
    fill(0, 0.25, 0.75)
    dot(pt, radius=1)

# y -= 25 from mesh
# cubic bez start, control control, end
warp = (
    Point(x=1, y=0),
    Point(x=43, y=-21),
    Point(x=85, y=13),
    Point(x=127, y=-8),
)
 
def flagWarpVec(pt, accuracy=0.001):
    # only active between warp[0].x and warp[-1].x
    if pt.x < warp[0].x or pt.x > warp[-1].x:
        print("Out of bounds:", pt)
        return pt, Point(0, 0)
    # Draw a line at pt.x through the warp; the pt on the warp is be the desired offset
    segments = splitCubic(*warp, pt.x, False)
    # no kinks/overlaps; there should be <=2 curves when split
    # we are really just interested in there being only one unique start/end at pt.x
    assert len(segments) <= 2, f"Too many segments: {segments}"
    wpt = Point(*(segments[0][-1]))
    deriv_pt = cubic_deriv_pos(1, *(Point(*s) for s in segments[0]))
    return Point(0, wpt.y), deriv_pt

def flagWarpPt(pt, accuracy=0.001):
    warpvec, derivvec = flagWarpVec(pt, accuracy)
    return (
        Point(pt.x + warpvec.x, pt.y + warpvec.y),
        Point(pt.x + derivvec.x, pt.y + derivvec.y),
    )
        

translate(y=-20)

stroke(None)
fill(0, 0, 0)
text("Shifts & Derivatives", (0, 6))

stroke(0.6, 0.6, 0.6)
line(warp[0], Point(warp[-1].x, warp[0].y))

for pt in warp:
    stroke(None)
    fill(0, 0.5, 0.75)
    dot(pt)

for x in range(1, 127 + 1, 4):
    pt = Point(x, 0)
    warp_pt, warp_deriv_vec = flagWarpPt(pt)

    stroke(0, 0.25, 0.9)
    strokeWidth(0.5)
    line(pt, warp_pt)

    stroke(0, 0.25, 0.9, 0.5)
    strokeWidth(0.25)
    if warp_deriv_vec.magnitude() != 0:
        tan1 = warp_deriv_vec.of_magnitude(1.5)
        tan2 = Point(-tan1.x, -tan1.y)
        line(tan1.add(warp_pt), tan2.add(warp_pt))

# OK, let's try to actually warp something?
translate(x=132, y=70)

stroke(None)
fill(0, 0, 0)
text("Original Shape(s)", (0, 64))

# we'll assume magic converted a rect into naive cubics
rect = (
    # bottom
    (
          Point(  0,  0),
          Point( 43,  0),
          Point( 85,  0),
          Point(128,  0),
    ),
    # rhs
    (
          Point(128,  0),
          Point(128, 20),
          Point(128, 40),
          Point(128, 60),
    ),
    # top
    (
          Point(128, 60),
          Point( 85, 60),
          Point( 43, 60),
          Point(  0, 60),
    ),
    # left
    (
          Point(  0, 60),
          Point(  0, 40),
          Point(  0, 20),
          Point(  0,  0),
    ),
)

fill(None)
stroke(1, 0, 0)
strokeWidth(1)
drawPath(bezOfCubics(rect))

for cubic in rect:
    for pt in cubic:
        stroke(None)
        fill(0, 0.5, 0.75)
        dot(pt)

translate(y=-20)

stroke(None)
fill(0, 0, 0)
text("Mangled Shape(s)", (2, 2))

translate(y=-64)

# original for reference
fill(None)
stroke(0.5, 0.5, 0.5)
strokeWidth(0.25)
drawPath(bezOfCubics(rect))

# warp

num_splits = 4
split_t = tuple(t / num_splits for t in range(1, num_splits))
split_mesh = []
for cubic in rect:
    for split_cubic in splitCubicAtT(*cubic, *split_t):
        split_mesh.append(tuple(Point(*pt) for pt in split_cubic))
split_mesh = tuple(split_mesh)

for cubic in split_mesh:
    for idx, pt in enumerate(cubic):
        stroke(None)
        fill(0, 0.5, 0.25 + 0.75)
        dot(pt, radius=1.5)

# mesh and warp modified slightly to cover entire [0, 128] on x
warp = (
    Point(x=0, y=0),
    Point(x=43, y=-21),
    Point(x=85, y=13),
    Point(x=128, y=-8),
)
mesh = (
  Point(  0, top+C),
  Point( 43, top-B+C),
  Point( 85, top+B-C),
  Point(128, top-C),
  Point(128, bot-C),
  Point( 85, bot+B-C),
  Point( 43, bot-B+C),
  Point(  0, bot+C),
)

warp_mesh = []
for cubic in split_mesh:
    # naive; just move according to shift
    # in the end expecting control point movement to differ from start/end point
    cubic = tuple(flagWarpPt(pt)[0] for pt in cubic)
    warp_mesh.append(cubic)
warp_mesh = tuple(warp_mesh)

fill(None)
stroke(0.75, 0.2, 0.2)
strokeWidth(0.5)
drawPath(bezOfCubics(warp_mesh))

for cubic in warp_mesh:    
    for idx, pt in enumerate(cubic):
        stroke(None)
        fill(0, 0.5, 0.25 + 0.25 * idx)
        dot(pt, radius=1.5)


