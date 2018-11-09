''' For converting a collection of hyperbolic right angled hexagons to a tilling of H^2. '''

from Queue import Queue
from surfaces import Tile
from hyperbolic import FLOAT_ERROR, dist_h, geodesic_centre, walk_along_geodesic, to_right_of, Isometry, interpolate
from hyperbolic import sqrt, fabs
import svg

next_side = [1,2,3,4,5,0]

class Tiling:
    def __init__(self, surface, radius=None, number=None, depth=None, S1=(-0.25, sqrt(3)/4), S2=(0.45, sqrt(3)/4)):
        ''' A tiling of H^2 from the hexagon decomposition of the surface. We start at S1 and go in
        the direction of S2 to place the first tile and tile out using tiles that;
         1) have a vertex within hyperbolic distance 'radius' from the centre,
         2) can be reached from the starting tile by passing through at most 'depth' number of tiles,
         3) are numbered less than 'number'.
        '''
        
        if radius is None and number is None and depth is None:
            raise ValueError('Cannot tile all of H^2')
        
        # Get the needed hexagons.
        self.hexagons = surface.hexagon_decomposition()
        
        self.tiles = [self.hexagons[0].basic_tile.move(0, 1, S1, S2)]
        tiles_vertices = set(self.tiles[0].approximate_vertices)
        
        Q = Queue()
        Q.put((self.tiles[0], 0))
        
        while not Q.empty():
            old_tile, old_depth = Q.get()
            
            # Skip tiles if the depth will be too far.
            if depth is not None and not old_depth < depth: continue
            
            old_hexagon = old_tile.parent
            basic_isometry = old_hexagon.basic_tile.get_move_isometry(0, 1, old_tile.vertices[0], old_tile.vertices[1])
            
            for i in range(6):
                # Skip unglued sides.
                if old_hexagon.gluing[i] is None: continue
                
                # Create the next new tile.
                new_tile = old_hexagon.adjacent_tile[i].move_isometry(basic_isometry)
                
                # Skip if the vertices of the tile are all too far from the centre.
                if radius is not None and not new_tile.min_distance() < radius: continue
                
                # Skip tiles we've already seen.
                if new_tile.approximate_vertices in tiles_vertices: continue
                
                self.tiles.append(new_tile)
                tiles_vertices.add(new_tile.approximate_vertices)
                Q.put((new_tile, old_depth+1))
                
                # Skip if we've created enough tiles. But we know that this test can now never be passed, so just give up.
                if number is not None and not len(self.tiles) < number: return
    
    def length_spectrum(self, base_lift=None, threshold=None):
        ''' Take a hexagonal tiling of a portion of H^2 and returns the
        minimum translation distance of any tile of the same type. If no
        base_lift tile is specified then the first one is used. A threshold
        may also be specified such that if a length is found with length
        less than it then the spectrum found so far is returned immediately.'''
        
        if base_lift is None: base_lift = self.tiles[0]
        base = base_lift.parent
        lengths = []
        for tile in self.tiles:
            if tile == base_lift: continue
            if tile.parent == base:
                h = Isometry(base_lift.vertices[0], base_lift.vertices[1], tile.vertices[0], tile.vertices[1])
                if h.type() == 2:  # Hyperbolic isometries only.
                    lengths.append(h.translation_length())
                    if threshold is not None and lengths[-1] < threshold-FLOAT_ERROR: return lengths
        
        return lengths
    
    def tile_containing(self, point):
        ''' Returns the tile containing point or None if no tile contains that point. '''
        
        for tile in self.tiles:
            if point in tile:
                return tile
        
        return None
    
    def orbit(self, point):
        ''' Returns a list of all the equivalent points in the tiling. '''
        
        base_tile = self.tile_containing(point)
        if base_tile is None: return None
        
        base = base_lift.parent
        
        orbit = []
        for tile in self.tiles:
            if tile.parent == base:
                h = Isometry(base_lift.vertices[0], base_lift.vertices[1], tile.vertices[0], tile.vertices[1])
                orbit.append(h(point))
        
        return orbit
    
    def epsilon_short_curve(self, epsilon):
        ''' Returns if the tiling contains a simple closed curve of length less than epsilon. '''
        
        return self.length_spectrum(threshold=epsilon)[-1] < epsilon-FLOAT_ERROR
    
    def systole_length(self):
        ''' Returns the length of the systole of the surface. '''
        
        return min(self.length_spectrum())
    
    def save_diagram(self, filename, size=1000):
        ''' Converts a collection of tiles to a scene of dimension size x size and saves it as a .svg file. '''
        s = size/2
        scene = svg.Scene('tiling', size, size)
        scene.add(svg.Circle((s,s),s,(255,255,255)))  # Boundary at infinity.
        
        # Draw each tile.
        for tile in self.tiles: scene.extend(tile.diagram(size))
        
        scene.write_svg(filename)

