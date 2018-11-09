import cmath
from math import cosh, sinh, tanh, acosh, sqrt, fabs
# from mpmath import cosh, sinh, acosh, sqrt, fabs  # Use the mpmath module to work to arbitary precision.

# Define some constants once.
infty = float('inf')
FLOAT_ERROR = 0.00005
boundary_circle = ((0.0,0.0), 1.0)

# Determine if two points are the same (within FLOAT_ERROR).
same_number = lambda m, n: -FLOAT_ERROR < (m - n) < FLOAT_ERROR
same_point = lambda m, n: same_number(m[0], n[0]) and same_number(m[1], n[1])

norm_2 = lambda u: u[0]*u[0] + u[1]*u[1]
norm = lambda u: sqrt(u[0]*u[0] + u[1]*u[1])
dist_h  = lambda u, v: acosh(1 + 2 * norm_2((u[0] - v[0], u[1] - v[1])) / ((1 - norm_2(u)) * (1 - norm_2(v))))

point_to_complex = lambda P: P[0]+1j*P[1]

class Isometry:
    def __init__(self, A, B, A2, B2):
        ''' Creates an isometry that sends A to A2 and B to a point on the geodesic [A2,B2].'''
        
        Ac, Bc, A2c, B2c = [point_to_complex(p) for p in [A,B,A2,B2]]
        
        t = cmath.sqrt((B2c-A2c)*(1-Bc*Ac.conjugate()) / ((Bc-Ac) * (1-B2c*A2c.conjugate())))
        # Later, rescaling will account for the fact that abs(t) != 1.
        a = t - t.conjugate() * Ac.conjugate() * A2c
        b = t.conjugate() * A2c - t * Ac
        
        # Rescale so the determinant is 1 and we have a matrix in PSL(2, C).
        d = sqrt(abs(a)**2 - abs(b)**2)
        a /= d
        b /= d
        
        # self.matrix = (a, b, b.conjugate(), a.conjugate())
        self.compressed_matrix = (a.real, a.imag, b.real, b.imag)
        self.trace = 2 * a.real  # Store the trace as well, tr() = a + a.conjugate() == 2*a.real
    
    def __call__(self, P):
        ''' Computes the action of the isometry on the point P.'''
        A =  self.compressed_matrix[2] + self.compressed_matrix[0]*P[0] - self.compressed_matrix[1]*P[1]
        B =  self.compressed_matrix[3] + self.compressed_matrix[1]*P[0] + self.compressed_matrix[0]*P[1]
        C =  self.compressed_matrix[0] + self.compressed_matrix[2]*P[0] + self.compressed_matrix[3]*P[1]
        D = -self.compressed_matrix[1] - self.compressed_matrix[3]*P[0] + self.compressed_matrix[2]*P[1]
        DD = C*C + D*D
        return ((A*C + B*D)/DD, (B*C - A*D)/DD)
    
    def type(self):
        ''' Returns the isometry type, 0 == elliptic, 1 == parabolic, 2 == hyperbolic.
        Note that NO checks are made that T is in the standard form (det == 1).'''
        
        trace_squared = self.trace**2
        
        if trace_squared < 4 - FLOAT_ERROR:  return 0
        if same_number(trace_squared, 4):    return 1
        if trace_squared > 4 + FLOAT_ERROR:  return 2
    
    def translation_length(self):
        ''' Returns the translation length of the hyperbolic isometry. Note that
        NO checks are made that T is hyperbolic and in the standard format (det == 1).'''
        
        return 2*acosh(abs(self.trace) / 2)

def vector(u, v):
    ''' Given u, v (Eucl) returns the vector from u to v. '''
    
    return (v[0]-u[0], v[1]-u[1])

def midpoint(u, v):
    ''' Given u, v (eucl) returns the midpoint of the line joining them. '''
    
    return ((u[0]+v[0])/2, (u[1]+v[1])/2)

def normalise(u):
    ''' Given u (Eucl) a none-zero vector returns the unit vector in that direction. '''
    
    norm_u = norm(u)
    return (u[0]/norm_u, u[1]/norm_u)

def form_positive_basis(U, V):
    ''' Give two (Eucl) vectors u & v, returns if u & v form a positive basis. '''
    
    return U[0] * V[1] >= U[1] * V[0]

def inside_circle(circle, p):
    ''' Given a circle (C, R) and a point p, returns if P is inside (C, R). '''
    (C, R) = circle
    
    return norm_2(vector(C, p)) < R * R

def to_right_of(p, u, v):
    ''' Returns if P is to the RHS of the geodesic connecting from u to v. '''
    
    C, R = geodesic_centre(u, v)
    U, V = vector(C, u), vector(C, v)
    
    return form_positive_basis(U, V) ^ inside_circle((C, R), p)

def perpendicular_vector(vector):
    ''' Given a vector in direction D starting at S, returns a vector perpendicuar 
    to D starting at S. '''
    (S, D) = vector
    
    return S, (D[1], -D[0])

def perpendicular_bisector(u, v):
    ''' Given u, v (eucl) returns the start and direction of the perpendicular bisector of [u,v]. '''
    
    return perpendicular_vector((midpoint(u, v), vector(u, v)))

def intersect_vectors(vector1, vector2):
    ''' Given 2 vectors, returns their point of intersection (or (infty,infty) if parallel).'''
    (S1,D1), (S2,D2) = vector1, vector2
    
    d = D1[0]*D2[1] - D1[1]*D2[0]
    if d == 0: return (infty, infty)
    
    A1 = (S2[0]*(S2[1]+D2[1]) - S2[1]*(S2[0]+D2[0]))
    A2 = (S1[0]*(S1[1]+D1[1]) - S1[1]*(S1[0]+D1[0]))
    
    Px = (D1[0]*A1-A2*D2[0]) / d
    Py = (D1[1]*A1-A2*D2[1]) / d

    return (Px, Py)

def intersect_circles(circle1, circle2):
    ''' Returns the two points of intersection of circles that meet.'''
    (C1, R1), (C2, R2) = circle1, circle2
    
    d = norm(vector(C1, C2))
    if d > R1 + R2 or d < fabs(R1-R2): return None
    if d == 0 and R1 == R2: return None
    
    a = (R1*R1 - R2*R2 + d*d) / (2*d)
    h = sqrt(R1*R1 - a*a)
    
    P2 = (C1[0] + (C2[0]-C1[0]) * a / d, C1[1] + (C2[1]-C1[1]) * a / d)
    
    if h == 0: return P2
    
    D = (h*(C2[1]-C1[1])/d, h*(C2[0]-C1[0])/d)
    
    return ((P2[0]+D[0], P2[1]-D[1]), (P2[0]-D[0],P2[1]+D[1]))

def involute_in_circle(P, circle=boundary_circle):
    ''' Returns the involution in (C,R) of P. '''
    (C,R) = circle
    
    Q = vector(C, P)
    t = norm_2(Q)
    return (Q[0]*R/t + C[0], Q[1]*R/t + C[1]) if t != 0 else (infty,infty)

def compute_circle(a, b, c):
    ''' Returns the (eulc) center and radius of the unique circle passing through A, B & C. '''
    # Should check for non-co-lineararity.
    
    C = intersect_vectors(perpendicular_bisector(a, b), perpendicular_bisector(b, c))
    R = norm(vector(C, a))
    
    # assert same_number(R, norm(vector(C, b))) and same_number(R, norm(vector(C, c)))
    
    return C, R

def boundary_geodesic_centre(u, v):
    ''' Given two points (eucl) on \partial H^2 returns the (eucl) centre & radius of the circle 
    passing through u and v perpendicular to \partial H^2. '''
    
    C = involute_in_circle(midpoint(u, v))
    R = norm((C[0]-u[0], C[1]-u[1]))
    
    return C, R

def geodesic_centre(u, v):
    ''' Given two points (eucl) in H^2 returns the (eucl) centre & radius of the circle 
    passing through u and v perpendicular to \partial H^2. '''
    
    norm_2_u, norm_2_v = norm_2(u), norm_2(v)
    
    if norm_2_u == 0 or norm_2_v == 0: return (infty,infty), infty
    if norm_2_u == 1 and norm_2_v == 1: return boundary_geodesic_centre(u, v)
    
    return compute_circle(u,v,involute_in_circle(u if norm_2(u) != 1 else v))

def geodesic_endpoints(u, v):
    ''' Given two points u and v returns the coordinates of the endpoints of the geodesic
    connecting them on \partial H^2. '''
    # !?! This should be MUCH smarter.
    
    geodesic = geodesic_centre(u,v)
    if geodesic == ((infty, infty), infty):
        return (normalise(u), normalise(v))
    else:
        return intersect_circles(geodesic, boundary_circle)

def perpendicular_circle(u, v):
    ''' Given 2 points u & v, computes the center & radius of the unique hyperbolic circle
    passsing through u, perpendicular to [u,v]. '''
    
    if norm_2(u) == 0 or norm_2(v) == 0: return (infty,infty), infty
    
    C1, R1 = geodesic_centre(u,v)
    
    C = intersect_vectors(perpendicular_bisector(u,involute_in_circle(u)), perpendicular_vector((u,vector(u, C1))))
    R = norm(vector(C, u))
    return C, R

def walk_along_geodesic(circle, u, t):
    ''' Returns the points on the circle (C,R) hyperbolic distance t from u, when u is a point on the circle. '''
    (C,R) = circle
    assert same_number(R, norm(vector(C, u)))  # Check u is actually on the geodesic.
    if same_number(t,0): return u, u
    
    return intersect_circles(hyperbolic_circle(u, t), (C,R))

def hyperbolic_circle(p, r):
    ''' Given a point p (eulcidean), and a hyperbolic radius r, returns the 
    Euclidean center and radius of the hyperbolic circle about p of hyp radius r. '''
    
    r = fabs(r)
    d = dist_h(p, (0,0))
    
    A = norm(p)
    t1, t2 = tanh(( r+d)/2) / A, tanh((-r+d)/2) / A
    
    t = (t1 + t2) / 2
    
    C = (p[0] * t, p[1] * t)
    R = fabs(t1 - t) * norm(p)
    return C, R

def interpolate(u, v, t):
    ''' Returns the point in H^2 that is a hyperbolic distance t along the geodesic [u,v]. '''
    uv_dist = dist_h(u, v)
    C1, R1 = hyperbolic_circle(u, t    )
    C2, R2 = hyperbolic_circle(v, uv_dist-t)
    
    v = vector(C2, C1)
    R = norm(v)
    
    return (C2[0] + R2/R * v[0], C2[1] + R2/R * v[1])

