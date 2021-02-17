from typing import NamedTuple

# from waveflag.c
top = 21
bot = 128-top
B = 21
C = 4

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


mesh = (
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

size(128, 196)

def dot(at, radius=2):
    oval(at.x-radius/2, at.y-radius/2, radius, radius)

# draw the basic shape and control points
#	cairo_line_to(cr,   M(0));
#	cairo_curve_to(cr,  M(1), M(2), M(3));
#	cairo_line_to(cr,   M(4));
#	cairo_curve_to(cr,  M(5), M(6), M(7));
#	cairo_close_path (cr);
path = BezierPath()
path.moveTo(mesh[0])
path.curveTo(mesh[1], mesh[2], mesh[3])
path.lineTo(mesh[4])
path.curveTo(mesh[5], mesh[6], mesh[7])
path.closePath()

fill(None)
stroke(1, 0, 0)
drawPath(path)

stroke(0.8, 0.8, 0.8)
strokeWidth(0.4)
line(mesh[0], mesh[1])
line(mesh[1], mesh[2])
line(mesh[2], mesh[3])

line(mesh[4], mesh[5])
line(mesh[5], mesh[6])
line(mesh[6], mesh[7])

# Oncurve points
fill(0, 0, 0)
stroke(None)
for idx in oncurve:
    dot(mesh[idx])

# Offcurve points
fill(0.5, 0.5, 0.5)
stroke(None)
for idx in offcurve:
    pt = mesh[idx]
    oval(pt.x-1, pt.y-1, 2, 2)

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
    

# points along curve

mesh_deriv = cubic_deriv(mesh[0], mesh[1], mesh[2], mesh[3])
print("cubic", mesh[0], mesh[1], mesh[2], mesh[3])
print("deriv", mesh_deriv)
line(mesh_deriv[0], mesh_deriv[1])
line(mesh_deriv[1], mesh_deriv[2])

stroke(None)
num_pts = 20
for t in range(0, num_pts + 1):
    pt = cubic_pos(t/num_pts, mesh[0], mesh[1], mesh[2], mesh[3])
    fill(0, 0, 0.75)
    dot(pt, radius=1.5)
    
    pt = quad_pos(t, *mesh_deriv)
    fill(0, 0.25, 0.75)
    dot(pt, radius=1)

# to "warp" a given point x,y based on insider knowledge of our specific cubic
# if x < min x, 0. if x > max x, 0
# find the nearest point on bezier
# on-curve movement is the difference between nearest.y and mesh[0].y
# warp bezier can place p0.y at 0 so on-curve move is just nearest.y

# y -= 25 from mesh
# cubic bez start, control control, end
warp = (
    Point(x=1, y=0),
    Point(x=43, y=-21),
    Point(x=85, y=13),
    Point(x=127, y=-8),
)
 
from fontTools.misc.bezierTools import splitCubic

# Take advantage of what we know about our warp:
# y-only
# only active between warp[0].x and warp[-1].x
# no kinks/overlaps; there should be only 2 curves when split
# TODO cache based on x; likely drawings reuse same x
def flagWarpVec(pt, accuracy=0.001):
    if pt.x < warp[0].x or pt.x > warp[-1].x:
        return pt, Point(0, 0)
    segments = splitCubic(*warp, pt.x, False)
    # lazy; really we care there is only one start/end at pt.x
    assert len(segments) == 2
    wpt = Point(*(segments[0][-1]))
    deriv_pt = cubic_deriv_pos(1, *(Point(*s) for s in segments[0]))
    return Point(0, wpt.y), deriv_pt

def flagWarpPt(pt, accuracy=0.001):
    warpvec, derivvec = flagWarpVec(pt, accuracy)
    return (
        Point(pt.x + warpvec.x, pt.y + warpvec.y),
        Point(pt.x + derivvec.x, pt.y + derivvec.y),
    )
        

translate(y=156)

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

    