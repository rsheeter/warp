"""
Python transcription of https://github.com/erich666/GraphicsGems/blob/master/gems/FitCurves.c
as it existed 2/24/2020.
"""
from math import sqrt
from typing import NamedTuple, Sequence, Tuple


class Point(NamedTuple):
    x: float
    y: float

    def magnitude(self) -> float:
        return sqrt(self.x * self.x + self.y * self.y)

    # V2AddII; it's weird this works on vectors...
    def add(self, other) -> "Point":
        return Point(self.x + other.x, self.y + other.y)

    # V2SubII
    def sub(self, other) -> "Point":
        return Point(self.x - other.x, self.y - other.y)

    def scale(self, s: float) -> "Point":
        return Point(self.x * s, self.y * s)

    @staticmethod
    def tan_left(points: Sequence["Point"]) -> "Vector":
        assert len(points) > 1
        return Vector.p0_to_p1(points[0], points[1]).normalize()

    @staticmethod
    def tan_right(points: Sequence["Point"]) -> "Vector":
        assert len(points) > 1
        return Vector.p0_to_p1(points[-1], points[-2]).normalize()

    @staticmethod
    def tan_center(points: Sequence["Point"]) -> "Vector":
        assert len(points) == 3
        v1 = Vector.p0_to_p1(points[1], points[0])
        v2 = Vector.p0_to_p1(points[2], points[1])
        tHatCenter = Vector((v1.x + v2.x) / 2, (v1.y + v2.y) / 2)
        return tHatCenter.normalize()


CubicBezier = Tuple[Point, Point, Point, Point]


class Vector(Point):
    def squared_length(self) -> float:
        return self.x * self.x + self.y * self.y

    def length(self) -> float:
        return sqrt(self.squared_length())

    def normalize(self) -> "Vector":
        length = self.length()
        x, y = 0.0, 0.0
        if length != 0.0:
            x = self.x / length
            y = self.y / length

        return Vector(x, y)

    def negate(self) -> "Vector":
        return Vector(-self.x, -self.y)

    @staticmethod
    def p0_to_p1(p0: Point, p1: Point) -> "Vector":
        return Vector(*(p1.sub(p0)))

    @staticmethod
    def dot_product(v0: "Vector", v1: "Vector") -> float:
        return v0.x * v1.x + v0.y * v1.y


def degree(curve):
    deg = len(curve) - 1
    assert deg in range(1, 4), f"we don't expect to handle degrees outside [1,3]: {deg}"
    return deg


def derivative(curve):
    curve_degree = degree(curve)
    deriv = []
    for i in range(curve_degree):
        deriv.append(
            Point(
                (curve[i + 1].x - curve[i].x) * curve_degree,
                (curve[i + 1].y - curve[i].y) * curve_degree,
            )
        )
    return tuple(deriv)


#  Reparameterize:
#     Given set of points and their parameterization, try to find
#     a better parameterization.
def _reparameterize(points, first, last, u, curve):
    new_ts = []
    for i in range(first, last + 1):
        new_ts.append(newton_raphson_root(curve, points[i], u[i - first]))
    return tuple(new_ts)


def newton_raphson_root(curve, pt, u):
    # Compute Q(u)
    curve_pt = point_at_t(curve, u)

    # Generate control vertices for curve'
    # TODO == cubic derivative; should just have derivative(curve)?
    curve_prime = derivative(curve)

    # Generate control vertices for curve''
    curve_prime_prime = derivative(curve_prime)

    # Compute Q'(u) and Q''(u)
    curve_prime_pt = point_at_t(curve_prime, u)
    curve_prime_prime_pt = point_at_t(curve_prime_prime, u)

    # Compute f(u)/f'(u)
    numerator = (curve_pt.x - pt.x) * (curve_prime_pt.x) + (curve_pt.y - pt.y) * (
        curve_prime_pt.y
    )
    denominator = (
        curve_prime_pt.x * curve_prime_pt.x
        + curve_prime_pt.y * curve_prime_pt.y
        + (curve_pt.x - pt.x) * curve_prime_prime_pt.x
        + (curve_pt.y - pt.y) * curve_prime_prime_pt.y
    )
    if denominator == 0.0:
        return u

    # u = u - f(u)/f'(u)
    u_prime = u - (numerator / denominator)
    return u_prime


# BezierII
def point_at_t(curve, t):
    curve_degree = degree(curve)
    mut_curve = list(curve)

    # Triangle computation
    for i in range(1, curve_degree + 1):
        for j in range(0, curve_degree - i + 1):
            mut_curve[j] = Point(
                (1.0 - t) * mut_curve[j].x + t * mut_curve[j + 1].x,
                (1.0 - t) * mut_curve[j].y + t * mut_curve[j + 1].y,
            )
    return mut_curve[0]


# B0, B1, B2, B3 :
# Bezier multipliers


def _B0(u: float) -> float:
    tmp = 1.0 - u
    return tmp * tmp * tmp


def _B1(u: float) -> float:
    tmp = 1.0 - u
    return 3 * u * (tmp * tmp)


def _B2(u: float) -> float:
    tmp = 1.0 - u
    return 3 * u * u * tmp


def _B3(u: float) -> float:
    return u * u * u


# GenerateBezier
def _generate_bezier(points, first, last, u_prime, tan_left, tan_right) -> CubicBezier:
    num_pts = last - first + 1  # number of pts in sub-curve

    # Precomputed rhs for eqn; Compute the A's
    A = tuple(
        (tan_left.scale(_B1(u_prime[i])), tan_right.scale(_B2(u_prime[i])))
        for i in range(num_pts)
    )

    # Create C and X matrices
    C = [[0.0, 0.0], [0.0, 0.0]]
    X = [0.0, 0.0]

    for i in range(num_pts):
        C[0][0] += Vector.dot_product(A[i][0], A[i][0])
        C[0][1] += Vector.dot_product(A[i][0], A[i][1])

        C[1][0] = C[0][1]
        C[1][1] += Vector.dot_product(A[i][1], A[i][1])

        vec_tmp = points[first + i].sub(
            points[first]
            .scale(_B0(u_prime[i]))
            .add(
                points[first]
                .scale(_B1(u_prime[i]))
                .add(
                    points[last]
                    .scale(_B2(u_prime[i]))
                    .add(points[last].scale(_B3(u_prime[i])))
                )
            )
        )

        X[0] += Vector.dot_product(A[i][0], vec_tmp)
        X[1] += Vector.dot_product(A[i][1], vec_tmp)

    # Compute the determinants of C and X
    det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1]
    det_C0_X = C[0][0] * X[1] - C[1][0] * X[0]
    det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1]

    # Finally, derive alpha values
    alpha_l = alpha_r = 0.0
    if det_C0_C1 != 0.0:
        alpha_l = det_X_C1 / det_C0_C1
        alpha_r = det_C0_X / det_C0_C1

    # If alpha negative, use the Wu/Barsky heuristic (see text)
    # (if alpha is 0, you get coincident control points that lead to
    # divide by zero in any subsequent NewtonRaphsonRootFind() call.
    seg_length = Vector.p0_to_p1(points[first], points[last]).length()
    epsilon = 1.0e-6 * seg_length
    if alpha_l < epsilon or alpha_r < epsilon:
        # fall back on standard (probably inaccurate) formula, and subdivide further if needed.
        dist = seg_length / 3.0
        return (
            points[first],
            points[first].add(tan_left.scale(dist)),
            points[last].add(tan_right.scale(dist)),
            points[last],
        )

    #  First and last control points of the Bezier curve are
    #  positioned exactly at the first and last data points
    #  Control points 1 and 2 are positioned an alpha distance out
    #  on the tangent vectors, left and right, respectively
    return (
        points[first],
        points[first].add(tan_left.scale(alpha_l)),
        points[last].add(tan_right.scale(alpha_r)),
        points[last],
    )


# ChordLengthParameterize
#    Assign parameter values to digitized points
#    using relative distances between points.
def _chord_length_parameterize(points, first, last) -> Sequence[float]:
    u = [0.0] * (last - first + 1)
    for i in range(first + 1, last + 1):
        u[i - first] = (
            u[i - first - 1] + Vector.p0_to_p1(points[i - 1], points[i]).length()
        )

    for i in range(first + 1, last + 1):
        u[i - first] = u[i - first] / u[last - first]

    return u


# TODO better "u" handling
def _max_error(points, first, last, curve, u) -> Tuple[float, int]:
    assert degree(curve) == 3
    max_err = 0.0
    split_at = int((last - first + 1) / 2)
    for i in range(first + 1, last):
        pt = point_at_t(curve, u[i - first])
        dist = Vector.p0_to_p1(points[i], pt).squared_length()
        if dist > max_err:
            max_err = dist
            split_at = i

    return max_err, split_at


# Fit a Bezier curve to a (sub)set of digitized points
def _fit_cubics(
    points, first, last, tan_left, tan_right, max_squared_error
) -> Tuple[CubicBezier]:
    iteration_error = max_squared_error * 4.0
    max_iterations = 4  # Max times to try iterating
    num_pts = last - first + 1

    assert num_pts > 1

    # Use heruistic for degenerate region
    if num_pts == 2:
        dist = Vector.p0_to_p1(points[0], points[1]).length() / 3.0
        return (
            (
                points[0],
                points[0].add(tan_left.scale(dist)),
                points[1].add(tan_right.scale(dist)),
                points[1],
            ),
        )

    #  Parameterize points, and attempt to fit curve
    u = _chord_length_parameterize(points, first, last)
    curve = _generate_bezier(points, first, last, u, tan_left, tan_right)

    # Find max deviation of points to fitted curve
    max_err, split_pt = _max_error(points, first, last, curve, u)
    if max_err < max_squared_error:
        return (curve,)

    # If error not too large, try some reparameterization
    # and iteration
    if max_err < iteration_error:
        for i in range(max_iterations):
            u_prime = _reparameterize(points, first, last, u, curve)
            curve = _generate_bezier(points, first, last, u_prime, tan_left, tan_right)
            max_err, split_pt = _max_error(points, first, last, curve, u_prime)
            if max_err < max_squared_error:
                return (curve,)
        u = u_prime

    # Fitting failed -- split at max error point and fit recursively
    tan_center = Point.tan_center(points[split_pt - 1 : split_pt + 2])
    return (
        *_fit_cubics(points, first, split_pt, tan_left, tan_center, max_err),
        *_fit_cubics(points, split_pt, last, tan_center.negate(), tan_right, max_err),
    )


def _as_point(pt):
    if isinstance(pt, Point):
        return pt
    if type(pt) == tuple and len(pt) == 2:
        return Point(*pt)
    raise ValueError(f"{type(pt)} does not convert to Point")


def fit_cubics(points, max_squared_error) -> Tuple[CubicBezier]:
    assert len(points) > 2, "No fun with < 3 points"
    points = tuple(_as_point(pt) for pt in points)
    tan_left = Point.tan_left(points)
    tan_right = Point.tan_right(points)
    return _fit_cubics(
        points, 0, len(points) - 1, tan_left, tan_right, max_squared_error
    )


def main():
    points = (
        Point(0.0, 0.0),
        Point(0.0, 0.5),
        Point(1.1, 1.4),
        Point(2.1, 1.6),
        Point(3.2, 1.1),
        Point(4.0, 0.2),
        Point(4.0, 0.0),
    )
    max_sqr_err = 4.0
    for cubic in fit_cubics(points, max_sqr_err):
        print("cubic")
        for pt in cubic:
            print(",".join("%.2f" % v for v in pt))


if __name__ == "__main__":
    main()
