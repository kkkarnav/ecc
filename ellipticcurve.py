class EllipticCurve(object):
    def __init__(self, a, b):
        # assume we're already in the Weierstrass form
        self.a = a
        self.b = b

        # if the discriminant is 0 then the curve is not smooth and thus invalid
        self.discriminant = -16 * (4 * a * a * a + 27 * b * b)
        if not self.discriminant:
            raise Exception(f"The curve {self} is not smooth!")

    # checks if a point is on the curve
    def test_point(self, x, y):
        return y * y == (x * x * x) + (self.a * x) + self.b

    # checks if two points on the curve are the same point
    def __eq__(self, other):
        return (self.a, self.b) == (other.a, other.b)

    def __str__(self):
        return f"y^2 = x^3 + {str(self.a)} + {str(self.b)}"


class Point(object):
    def __init__(self, curve, x, y):
        self.curve = curve  # the curve containing this point
        self.x = x
        self.y = y

        if not curve.test_point(x, y):
            raise Exception(f"The point {(self.x, self.y)} is not on the given curve {curve}")

    # reflects point about y
    def __neg__(self):
        return Point(self.curve, self.x, -self.y)

    # time complexity: O(n), space complexity: O(1) where n is
    def __add__(self, other):

        # point O at infinity will always sum up to itself
        if isinstance(other, Ideal):
            return self

        # if the two points are the same, take a tangent from the point which will intersect the curve at some point
        # x3 y3
        if (self.x, self.y) == (other.x, other.y):
            # if the points are distinct but y is 0, the intersection to the curve
            # will be at infinity and x3 y3 will be the ideal point
            if not self.y:
                return Ideal(self.curve, self.x, self.y)

            # slope of the tangent line
            slope_m = (3 * self.x * self.x + self.curve.a) / (2 * self.y)

        else:
            # if the points are distinct but the line intersecting them is vertical, the intersection to the curve
            # will be at infinity and x3 y3 will be the ideal point
            if self.x == other.x:
                return Ideal(self.curve, self.x, self.y)

            # slope of the line joining the two points extending upto the curve
            slope_m = (other.y - self.y) / (other.x - self.x)

        # the point of intersection of the extended line with the elliptic curve
        x3 = (slope_m * slope_m) - other.x - self.x
        y3 = slope_m * (x3 - self.x) + self.y

        # the point is reflected about the y-axis for summation
        return Point(self.curve, x3, -y3)

    # syntactic sugar
    def __sub__(self, Q):
        return self + -Q

    # efficient scaling of the point (using an integer, not a field element)
    def __mul__(self, n):

        if not isinstance(n, int):
            raise Exception("Can't scale a point by something which isn't an int!")

        else:
            if n < 0:
                return -self * -n
            # evaluates to a point at the end of the curve, which is O
            if n == 0:
                return Ideal(self.curve, self.x, self.y)
            else:
                P = self
                R = self if n % 2 == 1 else Ideal(self.curve, self.x, self.y)

                # a power of two, allowing us to represent x as b_k for k from 0 until 2^k = i > n
                i = 2
                while i <= n:
                    # adding the point to itself
                    P += P

                    # add P to the product if its i^th bit is 1 (since that indicates that 1 x 2^i is part of nQ)
                    if n & i == i:
                        R += P

                    # multiply i by 2 by bitshifting in the binary representation
                    i <<= 1

                return R

    # syntactic sugar
    def __rmul__(self, n):
        return self * n

    # resolves weird issue where normal equality operator wasn't saying two points on the same curve with same (x,
    # y) weren't being evaluated as equal
    def __eq__(self, other):
        return (self.x, self.y, self.curve) == (other.x, other.y, other.curve)

    def __list__(self):
        return [self.x, self.y]


# The point O at the end of the curve
class Ideal(Point):
    def __init__(self, curve, x, y):
        Point.__init__(self, curve, x, y)
        self.curve = curve

    def __str__(self):
        return "Ideal"

    def __neg__(self):
        return self

    def __add__(self, other):
        return other

    def __mul__(self, n):
        if not isinstance(n, int):
            raise Exception("Can't scale a point by something which isn't an int!")
        else:
            return self


if __name__ == "__main__":
    curve = (EllipticCurve(a=-2, b=4))
    point1 = (Point(curve=curve, x=3, y=5))
    point2 = Point(curve=curve, x=-2, y=0)
    point3 = point2 - point1
    print(point3.__list__())
