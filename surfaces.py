from hyperbolic import cosh, sinh, acosh, sqrt, infty
from hyperbolic import to_right_of, walk_along_geodesic, perpendicular_circle, geodesic_centre, Isometry, dist_h, same_point, interpolate
from itertools import combinations
import svg

next_side = [1,2,3,4,5,0]

def smart_mod(x):
    ''' Returns x % 1.0 but over the interval (-0.5,0.5]. '''
    
    c = x % 1.0
    return c - 1 if c > 0.5 else c

class Tile:
    def __init__(self, vertices, parent):
        self.vertices = vertices
        self.approximate_vertices = tuple([(int(x[0] * 1000), int(x[1] * 1000)) for x in self.vertices])
        self.parent = parent
    def vertex_pairs(self):
        for i in range(6):
            yield self.vertices[i-1], self.vertices[i]
    def diameter(self):
        return max(dist_h(u,v) for u, v in combinations(self.vertices, 2))
    def get_move_isometry(self, index1, index2, S1, S2):
        return Isometry(self.vertices[index1], self.vertices[index2], S1, S2)
    def move_isometry(self, h):
        return Tile(tuple(map(h, self.vertices)), self.parent)
    def move(self, index1, index2, S1, S2):
        return self.move_isometry(self.get_move_isometry(index1, index2, S1, S2))
    def min_distance(self, p=(0,0)):
        return min(dist_h(u, p) for u in self.vertices)
    def max_distance(self, p=(0,0)):
        return max(dist_h(u, p) for u in self.vertices)
    def __contains__(self, p):
        return all(to_right_of(p,v1,v2) for v1, v2 in self.vertex_pairs())
    def diagram(self, size=1000):
        s = size / 2
        for v1, v2 in self.vertex_pairs():
            C, R = geodesic_centre(v1,v2)
            yield svg.Arc((s+v1[0]*s,s+v1[1]*s),(s+v2[0]*s,s+v2[1]*s),R*s,(0,0,0),not to_right_of(C, v1, v2))

class Right_Hexagon:
    def __init__(self, *side_lengths):
        assert len(side_lengths) == 6
        self.side_lengths = side_lengths
        
        self.gluing = [None] * 6
        self.adjacent_tile = [None] * 6
        self.basic_tile = self.place_tile()
    
    def set_gluing(self, side, target, target_side, offset):
        self.gluing[side] = (target, target_side, offset)
        target.gluing[target_side] = (self, side, offset)
        
        a, b = self.basic_tile.vertices[side], self.basic_tile.vertices[next_side[side]]
        a, b = interpolate(a, b, offset*self.side_lengths[side]), interpolate(b, a, -offset*self.side_lengths[side])
        self.adjacent_tile[side] = target.basic_tile.move(target_side, next_side[target_side], b, a)
        
        a, b = target.basic_tile.vertices[target_side], target.basic_tile.vertices[next_side[target_side]]
        a, b = interpolate(a, b, offset*target.side_lengths[target_side]), interpolate(b, a, -offset*target.side_lengths[target_side])
        target.adjacent_tile[target_side] = self.basic_tile.move(side, next_side[side], b, a)
    
    def place_tile(self, S1=(-0.25, sqrt(3)/4), S2=(0.45, sqrt(3)/4), start=0):
        ''' Prduces a tile of this hexagon starting at S1, going in the direction of S2 and 
        turning right at each vertex, starts with side start (by default side 0). '''
        
        sides = self.side_lengths[start:] + self.side_lengths[:start]
        P = [None] * 7
        
        # We build P[0] and P[1] manually.
        P[0] = S1
        P1, P2 = walk_along_geodesic(geodesic_centre(S1, S2), P[0], sides[0])
        P[1] = P1 if dist_h(S2, P1) < dist_h(S2, P2) else P2
        
        for i in range(1,6):
            P1, P2 = walk_along_geodesic(perpendicular_circle(P[i], P[i-1]), P[i], sides[i])
            P[i+1] = P1 if to_right_of(P1, P[i-1], P[i]) else P2  # Work out which of P1 & P2 we want.
        
        assert same_point(P[0], P[6])  # Make sure that we've roughly closed up.
        
        return Tile(tuple(P[6-start:6]+P[0:6-start]), self)
    
    def diameter(self):
        return self.basic_tile.diameter()

# Each pair of pants decomposes into 2 right angled hexagons, here 
# is the front one. The back one is its vertical mirror image.
# 
#             cuff[0]/2
#             ---------
#            /         \
# seem[1]   /           \ seem[2]
#          /             \
#         /               \
#         \               /
#          \             /
# cuff[2]/2 \           / cuff[1]/2
#            \         /
#             ---------
#              seem[0]

class Pants:
    def __init__(self, *cuff_lengths):
        assert cuff_lengths is None or len(cuff_lengths) == 3
        self.cuff_lengths = [0.01, 0.01, 0.01] if cuff_lengths is None else list(cuff_lengths)
        self.gluing = [None] * 3  # Each is a dict of the form {target, target_cuff, offset} and starts unglued.
        self.front_hexagon = None
        self.back_hexagon = None
    
    def set_gluing(self, cuff, target, target_cuff, offset):
        self.gluing[cuff] = {'target':target, 'target_cuff':target_cuff, 'offset':offset}
        target.gluing[target_cuff] = {'target':self, 'target_cuff':cuff, 'offset':offset}
        target.cuff_lengths[target_cuff] = self.cuff_lengths[cuff]
    
    def set_cuff(self, cuff, new_cuff_length, new_cuff_offset):
        self.cuff_lengths[cuff] = new_cuff_length
        
        if self.gluing[cuff] != None:
            self.gluing[cuff]['offset'] = new_cuff_offset
            
            target = self.gluing[cuff]['target']
            target_cuff = self.gluing[cuff]['target_cuff']
            target.cuff_lengths[target_cuff] = new_cuff_length
            target.gluing[target_cuff]['offset'] = new_cuff_offset
    
    def compute_seam_lengths(self):
        c = map(lambda x: x/2, self.cuff_lengths)
        x = acosh( (cosh(c[1]) * cosh(c[2]) + cosh(c[0])) / (sinh(c[1]) * sinh(c[2])) ) if c[1] != 0 and c[2] != 0 else infty
        y = acosh( (cosh(c[2]) * cosh(c[0]) + cosh(c[1])) / (sinh(c[2]) * sinh(c[0])) ) if c[2] != 0 and c[0] != 0 else infty
        z = acosh( (cosh(c[0]) * cosh(c[1]) + cosh(c[2])) / (sinh(c[0]) * sinh(c[1])) ) if c[0] != 0 and c[1] != 0 else infty
        return x,y,z
    
    def hexagon_decompose(self):
        c = map(lambda x: x/2, self.cuff_lengths)
        s = self.compute_seam_lengths()
        
        # Make 2 new hexagons of the right size.
        front = Right_Hexagon(c[0], s[2], c[1], s[0], c[2], s[1])
        back = Right_Hexagon(c[0], s[1], c[2], s[0], c[1], s[2])
        
        # Glue these hexagons together.
        front.set_gluing(1, back, 5, 0.0)
        front.set_gluing(3, back, 3, 0.0)
        front.set_gluing(5, back, 1, 0.0)
        self.front_hexagon = front
        self.back_hexagon = back

class PantsDecomposition:
    def __init__(self, cuff_lengths, cuff_gluings=None):
        ''' Gluings are of the form (source_pant_index, source_pant_cuff, target_pant_index, target_pant_cuff, torsion).'''
        self.my_pants = [Pants(*cuff_length) for cuff_length in cuff_lengths]
        if cuff_gluings is not None: 
            for cuff_gluing in cuff_gluings:
                self.glue_cuffs(*cuff_gluing)
    
    def glue_cuffs(self, source_pant_index, source_pant_cuff, target_pant_index, target_pant_cuff, torsion):
        self.my_pants[source_pant_index].set_gluing(source_pant_cuff, self.my_pants[target_pant_index], target_pant_cuff, torsion)
    
    def unglue_cuffs(self, pant_index, pant_cuff):
        ''' Cuffs are of the form (pant_index, pant_cuff). '''
        for cuff in cuffs:
            if self.my_pants[pant_index].gluing[pant_cuff] is not None:
                target = self.my_pants[pant_index].gluing[pant_cuff]
                target['target'].gluing[target['target_cuff']] = None
            self.my_pants[pant_index].gluing[pant_cuff] = None
    
    def hexagon_decomposition(self):
        ''' Decomposes the surface into a collection of right angled hexagons, 2 for each pair of pants. 
        These can then be used to create a tiling of H^2 for example. '''
        for pant in self.my_pants: pant.hexagon_decompose()
        corresponding = {0:0, 2:4, 4:2}
        
        for pant in self.my_pants:
            for i in [0,1,2]:
                if pant.gluing[i] is not None:
                    pant.front_hexagon.set_gluing( i*2,
                        pant.gluing[i]['target'].front_hexagon, 
                        pant.gluing[i]['target_cuff']*2,
                        smart_mod(pant.gluing[i]['offset']) * 2)
                    pant.back_hexagon.set_gluing( corresponding[i*2],
                        pant.gluing[i]['target'].back_hexagon, 
                        corresponding[pant.gluing[i]['target_cuff']*2],
                        smart_mod(pant.gluing[i]['offset']) * 2)
        
        return [pant.front_hexagon for pant in self.my_pants] + [pant.back_hexagon for pant in self.my_pants]
    
    def diameter(self):
        return sum(h.diameter() for h in self.hexagon_decomposition())

    def information(self):
        boundary_components = len([1 for p in self.my_pants for c in p.gluing if c is None])
        Euler_characteristic = -len(self.my_pants)
        genus = (2 - Euler_characteristic - boundary_components) / 2
        return genus, boundary_components

