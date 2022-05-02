# credits Jeremy Kun, @j2kun

import fractions


# memoize calls to the class constructors for fields
# this helps typechecking by never creating two separate
# instances of a number class.
def memoize(f):
   cache = {}

   def memoizedFunction(*args, **kwargs):
      argTuple = args + tuple(kwargs)
      if argTuple not in cache:
         cache[argTuple] = f(*args, **kwargs)
      return cache[argTuple]

   memoizedFunction.cache = cache
   return memoizedFunction


# type check a binary operation, and silently typecast 0 or 1
def typecheck(f):
   def newF(self, other):
      if (hasattr(other.__class__, 'operatorPrecedence') and
            other.__class__.operatorPrecedence > self.__class__.operatorPrecedence):
         return NotImplemented

      if type(self) is not type(other):
         try:
            other = self.__class__(other)
         except TypeError:
            message = 'Not able to typecast %s of type %s to type %s in function %s'
            raise TypeError(message % (other, type(other).__name__, type(self).__name__, f.__name__))
         except Exception as e:
            message = 'Type error on arguments %r, %r for functon %s. Reason:%s'
            raise TypeError(message % (self, other, f.__name__, type(other).__name__, type(self).__name__, e))

      return f(self, other)

   return newF



# require a subclass to implement +-* neg and to perform typechecks on all of
# the binary operations finally, the __init__ must operate when given a single
# argument, provided that argument is the int zero or one
class DomainElement(object):
   operatorPrecedence = 1

   # the 'r'-operators are only used when typecasting ints
   def __radd__(self, other): return self + other
   def __rsub__(self, other): return -self + other
   def __rmul__(self, other): return self * other

   # square-and-multiply algorithm for fast exponentiation
   def __pow__(self, n):
      if type(n) is not int:
         raise TypeError

      Q = self
      R = self if n & 1 else self.__class__(1)

      i = 2
      while i <= n:
         Q = (Q * Q)

         if n & i == i:
            R = (Q * R)

         i = i << 1

      return R


   # requires the additional % operator (i.e. a Euclidean Domain)
   def powmod(self, n, modulus):
      if type(n) is not int:
         raise TypeError

      Q = self
      R = self if n & 1 else self.__class__(1)

      i = 2
      while i <= n:
         Q = (Q * Q) % modulus

         if n & i == i:
            R = (Q * R) % modulus

         i = i << 1

      return R



# additionally require inverse() on subclasses
class FieldElement(DomainElement):
   def __truediv__(self, other): return self * other.inverse()
   def __rtruediv__(self, other): return self.inverse() * other
   def __div__(self, other): return self.__truediv__(other)
   def __rdiv__(self, other): return self.__rtruediv__(other)

# strip all copies of elt from the end of the list
def strip(L, elt):
   if len(L) == 0: return L

   i = len(L) - 1
   while i >= 0 and L[i] == elt:
      i -= 1

   return L[:i+1]


# create a polynomial with coefficients in a field; coefficients are in
# increasing order of monomial degree so that, for example, [1,2,3]
# corresponds to 1 + 2x + 3x^2
@memoize
def polynomialsOver(field=fractions.Fraction):

   class Polynomial(DomainElement):
      operatorPrecedence = 2

      @classmethod
      def factory(cls, L):
         return Polynomial([cls.field(x) for x in L])

      def __init__(self, c):
         if type(c) is Polynomial:
            self.coefficients = c.coefficients
         elif isinstance(c, field):
            self.coefficients = [c]
         elif not hasattr(c, '__iter__') and not hasattr(c, 'iter'):
            self.coefficients = [field(c)]
         else:
            self.coefficients = c

         self.coefficients = strip(self.coefficients, field(0))
         self.indeterminate = 't'


      def isZero(self): return self.coefficients == []

      def __repr__(self):
         if self.isZero():
            return '0'

         return ' + '.join(['%s %s^%d' % (a, self.indeterminate, i) if i > 0 else '%s'%a
                              for i,a in enumerate(self.coefficients)])


      def __abs__(self): return len(self.coefficients) # the valuation only gives 0 to the zero polynomial, i.e. 1+degree
      def __len__(self): return len(self.coefficients)
      def __sub__(self, other): return self + (-other)
      def __iter__(self): return iter(self.coefficients)
      def __neg__(self): return Polynomial([-a for a in self])

      def iter(self): return self.__iter__()
      def leadingCoefficient(self): return self.coefficients[-1]
      def degree(self): return abs(self) - 1

      @typecheck
      def __eq__(self, other):
         return self.degree() == other.degree() and all([x==y for (x,y) in zip(self, other)])

      @typecheck
      def __ne__(self, other):
          return self.degree() != other.degree() or any([x!=y for (x,y) in zip(self, other)])

      @typecheck
      def __add__(self, other):
         newCoefficients = [sum(x) for x in zip_longest(self, other, fillvalue=self.field(0))]
         return Polynomial(newCoefficients)


      @typecheck
      def __mul__(self, other):
         if self.isZero() or other.isZero():
            return Zero()

         newCoeffs = [self.field(0) for _ in range(len(self) + len(other) - 1)]

         for i,a in enumerate(self):
            for j,b in enumerate(other):
               newCoeffs[i+j] += a*b

         return Polynomial(newCoeffs)


      @typecheck
      def __divmod__(self, divisor):
         quotient, remainder = Zero(), self
         divisorDeg = divisor.degree()
         divisorLC = divisor.leadingCoefficient()

         while remainder.degree() >= divisorDeg:
            monomialExponent = remainder.degree() - divisorDeg
            monomialZeros = [self.field(0) for _ in range(monomialExponent)]
            monomialDivisor = Polynomial(monomialZeros + [remainder.leadingCoefficient() / divisorLC])

            quotient += monomialDivisor
            remainder -= monomialDivisor * divisor

         return quotient, remainder


      @typecheck
      def __truediv__(self, divisor):
         if divisor.isZero():
            raise ZeroDivisionError
         return divmod(self, divisor)[0]


      @typecheck
      def __mod__(self, divisor):
         if divisor.isZero():
            raise ZeroDivisionError
         return divmod(self, divisor)[1]


   def Zero():
      return Polynomial([])


   Polynomial.field = field
   Polynomial.__name__ = '(%s)[x]' % field.__name__
   Polynomial.englishName = 'Polynomials in one variable over %s' % field.__name__
   return Polynomial

# so all IntegersModP are instances of the same base class
class _Modular(FieldElement):
   pass


@memoize
def IntegersModP(p):
   # assume p is prime

   class IntegerModP(_Modular):
      def __init__(self, n):
         try:
            self.n = int(n) % IntegerModP.p
         except:
            raise TypeError("Can't cast type %s to %s in __init__" % (type(n).__name__, type(self).__name__))

         self.field = IntegerModP

      @typecheck
      def __add__(self, other):
         return IntegerModP(self.n + other.n)

      @typecheck
      def __sub__(self, other):
         return IntegerModP(self.n - other.n)

      @typecheck
      def __mul__(self, other):
         return IntegerModP(self.n * other.n)

      def __neg__(self):
         return IntegerModP(-self.n)

      @typecheck
      def __eq__(self, other):
         return isinstance(other, IntegerModP) and self.n == other.n

      @typecheck
      def __ne__(self, other):
         return isinstance(other, IntegerModP) is False or self.n != other.n

      @typecheck
      def __divmod__(self, divisor):
         q,r = divmod(self.n, divisor.n)
         return (IntegerModP(q), IntegerModP(r))

      def inverse(self):
         # need to use the division algorithm *as integers* because we're
         # doing it on the modulus itself (which would otherwise be zero)
         x,y,d = extendedEuclideanAlgorithm(self.n, self.p)

         if d != 1:
            raise Exception("Error: p is not prime in %s!" % (self.__name__))

         return IntegerModP(x)

      def __abs__(self):
         return abs(self.n)

      def __str__(self):
         return str(self.n)

      def __repr__(self):
         return '%d (mod %d)' % (self.n, self.p)

      def __int__(self):
         return self.n

   IntegerModP.p = p
   IntegerModP.__name__ = 'Z/%d' % (p)
   IntegerModP.englishName = 'IntegersMod%d' % (p)
   return IntegerModP

# a general Euclidean algorithm for any number type with
# a divmod and a valuation abs() whose minimum value is zero
def gcd(a, b):
   if abs(a) < abs(b):
      return gcd(b, a)

   while abs(b) > 0:
      _,r = divmod(a,b)
      a,b = b,r

   return a


# extendedEuclideanAlgorithm: int, int -> int, int, int
# input (a,b) and output three numbers x,y,d such that ax + by = d = gcd(a,b).
# Works for any number type with a divmod and a valuation abs()
# whose minimum value is zero
def extendedEuclideanAlgorithm(a, b):
   if abs(b) > abs(a):
      (x,y,d) = extendedEuclideanAlgorithm(b, a)
      return (y,x,d)

   if abs(b) == 0:
      return (1, 0, a)

   x1, x2, y1, y2 = 0, 1, 1, 0
   while abs(b) > 0:
      q, r = divmod(a,b)
      x = x2 - q*x1
      y = y2 - q*y1
      a, b, x2, x1, y2, y1 = b, r, x1, x, y1, y

   return (x2, y2, a)



# isIrreducible: Polynomial, int -> bool
# determine if the given monic polynomial with coefficients in Z/p is
# irreducible over Z/p where p is the given integer
# Algorithm 4.69 in the Handbook of Applied Cryptography
def isIrreducible(polynomial, p):
   ZmodP = IntegersModP(p)
   if polynomial.field is not ZmodP:
      raise TypeError("Given a polynomial that's not over %s, but instead %r" %
                        (ZmodP.__name__, polynomial.field.__name__))

   poly = polynomialsOver(ZmodP).factory
   x = poly([0,1])
   powerTerm = x
   isUnit = lambda p: p.degree() == 0

   for _ in range(int(polynomial.degree() / 2)):
      powerTerm = powerTerm.powmod(p, polynomial)
      gcdOverZmodp = gcd(polynomial, powerTerm - x)
      if not isUnit(gcdOverZmodp):
         return False

   return True


# generateIrreduciblePolynomial: int, int -> Polynomial
# generate a random irreducible polynomial of a given degree over Z/p, where p
# is given by the integer 'modulus'. This algorithm is expected to terminate
# after 'degree' many irreducilibity tests. By Chernoff bounds the probability
# it deviates from this by very much is exponentially small.
def generateIrreduciblePolynomial(modulus, degree):
   Zp = IntegersModP(modulus)
   Polynomial = polynomialsOver(Zp)

   while True:
      coefficients = [Zp(random.randint(0, modulus-1)) for _ in range(degree)]
      randomMonicPolynomial = Polynomial(coefficients + [Zp(1)])
      print(randomMonicPolynomial)

      if isIrreducible(randomMonicPolynomial, modulus):
         return randomMonicPolynomial


# create a type constructor for the finite field of order p^m for p prime, m >= 1
@memoize
def FiniteField(p, m, polynomialModulus=None):
   Zp = IntegersModP(p)
   if m == 1:
      return Zp

   Polynomial = polynomialsOver(Zp)
   if polynomialModulus is None:
      polynomialModulus = generateIrreduciblePolynomial(modulus=p, degree=m)

   class Fq(FieldElement):
      fieldSize = int(p ** m)
      primeSubfield = Zp
      idealGenerator = polynomialModulus
      operatorPrecedence = 3

      def __init__(self, poly):
         if type(poly) is Fq:
            self.poly = poly.poly
         elif type(poly) is int or type(poly) is Zp:
            self.poly = Polynomial([Zp(poly)])
         elif isinstance(poly, Polynomial):
            self.poly = poly % polynomialModulus
         else:
            self.poly = Polynomial([Zp(x) for x in poly]) % polynomialModulus

         self.field = Fq

      @typecheck
      def __add__(self, other): return Fq(self.poly + other.poly)
      @typecheck
      def __sub__(self, other): return Fq(self.poly - other.poly)
      @typecheck
      def __mul__(self, other): return Fq(self.poly * other.poly)
      @typecheck
      def __eq__(self, other): return isinstance(other, Fq) and self.poly == other.poly

      def __pow__(self, n): return Fq(pow(self.poly, n))
      def __neg__(self): return Fq(-self.poly)
      def __abs__(self): return abs(self.poly)
      def __repr__(self): return repr(self.poly) + ' \u2208 ' + self.__class__.__name__

      @typecheck
      def __divmod__(self, divisor):
         q,r = divmod(self.poly, divisor.poly)
         return (Fq(q), Fq(r))


      def inverse(self):
         if self == Fq(0):
            raise ZeroDivisionError

         x,y,d = extendedEuclideanAlgorithm(self.poly, self.idealGenerator)
         if d.degree() != 0:
            raise Exception('Somehow, this element has no inverse! Maybe intialized with a non-prime?')

         return Fq(x) * Fq(d.coefficients[0].inverse())


   Fq.__name__ = 'F_{%d^%d}' % (p,m)
   return Fq
