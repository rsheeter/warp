from typing import NamedTuple

# from waveflag.c
top = 21
bot = 128-top
B = 21
C = 4

class Point(NamedTuple):
    x: float
    y: float

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

# Take advantage of what we know about our warp:
# y-only
# only active between warp[0].x and warp[-1].x
# no kinks/overlaps; safe to bsearch between min/max x
# TODO cache based on x; likely drawings reuse same x
def flagWarpVec(pt, accuracy=0.001):
    if pt.x < warp[0].x or pt.x > warp[-1].x:
        return pt, Point(0, 0)
    # binary search over t seeking same x as pt
    t = 1.0
    step = 0.5
    for i in range(100):
        wpt = cubic_pos(t, *warp)
        if abs(wpt.x - pt.x) <= accuracy:
            return (Point(0, wpt.y), cubic_deriv_pos(t, *warp))
        if wpt.x > pt.x:
            t -= step
        else:
            t += step
        assert t >= 0 and t <= 1
        step = step / 2
        
    raise ValueError(f"Unable to find warp for {pt}")

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
    warp_pt, warp_deriv_pt = flagWarpPt(pt)

    stroke(0, 0.25, 0.9)
    strokeWidth(0.5)
    line(pt, warp_pt)

    stroke(0, 0.25, 0.9, 0.5)
    strokeWidth(0.25)
    #line(pt, warp_deriv_pt)

translate(y=32)

for x in range(0, 120+1, 8):
    stroke(0.4, 0.4, 0.4)
    line((x, 0), (x+2, 0))
    
    warp_pt1, warp_deriv_pt1 = flagWarpPt(Point(x, 0))
    warp_pt2, warp_deriv_pt2 = flagWarpPt(Point(x + 2, 0))
    line((x + 0, 0), (x + 2, 0))
    line(warp_pt1, warp_pt2)

    stroke(0.7, 0.7, 0.7)
    warp_pt1, warp_deriv_pt1 = flagWarpPt(Point(x + 6, 0))
    warp_pt2, warp_deriv_pt2 = flagWarpPt(Point(x + 8, 0))
    line((x + 6, 0), (x+8, 0))
    line(warp_pt1, warp_pt2)
    