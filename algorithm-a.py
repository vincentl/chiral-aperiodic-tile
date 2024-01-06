import math
from functools import reduce
import argparse
from svgpathtools import Path, Line, QuadraticBezier, CubicBezier, Arc, wsvg

a = 1
b = 1
c = math.cos(math.pi/3)
s = math.sin(math.pi/3)

cmplx = lambda x,y: x + y*1j
zround = lambda c, p: round(64+c.real,p) + 1j*round(64+c.imag,p)

moves = [cmplx(0,0), 
         cmplx(c*b, s*b),
         cmplx(b, 0),
         cmplx(0,a),
         cmplx(s*a, c*a),
         cmplx(c*b, -s*b),
         cmplx(-c*b, -s*b),
         cmplx(s*a, -c*a),
         cmplx(0,-a),
         cmplx(0,-a),
         cmplx(-s*a, -c*a),
         cmplx(-c*b, s*b),
         cmplx(-b, 0),
         cmplx(0, a),
         cmplx(-s*a, c*a)]

# https://math.stackexchange.com/questions/3685824/a-mapping-that-converts-a-line-segment-to-another-one
# xA  = f[0].real
# yA  = f[0].imag
# xA' = f[1].real
# yA' = f[1].imag
# xB  = t[0].real
# yB  = t[0].imag
# xB' = t[1].real
# yB' = t[1].imag

def isometry(f,t):
    Ax = f[1].real - f[0].real
    Ay = f[1].imag - f[0].imag
    A2 = Ax*Ax + Ay*Ay
    Bx = t[1].real - t[0].real
    By = t[1].imag - t[0].imag
    ux = (Ax*Bx + Ay*By)/A2
    uy = (Ax*By - Ay*Bx)/A2
    tx = t[0].real - (f[0].real*(Ax*Bx + Ay*By) + f[0].imag*(Ay*Bx - Ax*By))/A2
    ty = t[0].imag - (f[0].real*(Ax*By - Ay*Bx) + f[0].imag*(Ax*Bx + Ay*By))/A2    
    return ([[ux, -uy],[uy, ux]], tx + 1j*ty)

def apply(I, pt):
    m, t = I
    x, y = pt.real, pt.imag
    x, y = m[0][0] * x + m[0][1] * y, m[1][0] * x + m[1][1] * y
    pt   = cmplx(x,y) + t
    return pt

def bezier_isometry(src, dst):
    return isometry([src.end, src.start], [dst.start, dst.end])

def bezier_rotate_isometry(edge, point, angle):
    rotate = math.cos(-angle) + 1j*math.sin(-angle) 
    src = [edge.start, edge.end]
    dst = [(p - point)*rotate + point for p in src]
    return isometry(src, dst)

def bezier_apply(I, curve):
    return CubicBezier(apply(I, curve.start), apply(I, curve.control1), apply(I, curve.control2), apply(I, curve.end))

# CubicBezier(1+0j, 0.8+0.5j, 0.1+0.2j, 0+0j)
# base curve defined from (0,0) to (1,0)
# want curve from s to e
def bezier_A(s,e):
    I = isometry([0j,1+0j],(s,e))
    c0 = apply(I, 0.1 + 0.2j)
    c1 = apply(I, 0.8 + 0.5j)
    #c0 = apply(I, 0.1 + 0.01j)
    #c1 = apply(I, 0.9 + 0.01j)
    return CubicBezier(s, c0, c1, e)

def bezier_B(s,e):
    I = isometry([0j,1+0j],(s,e))
    c0 = apply(I, 0.2 - 0.5j)
    c1 = apply(I, 0.9 - 0.2j)
    #c0 = apply(I, 0.1 - 0.01j)
    #c1 = apply(I, 0.9 - 0.01j)
    return CubicBezier(s, c0, c1, e)

def bezier(s,e,which):
    f = [bezier_B,bezier_A][which % 2]
    return f(s,e)

marks = lambda i: (f'{zround(i.start,3):.3f}:{zround(i.end,3):.3f}', f'{zround(i.end,3):.3f}:{zround(i.start,3):.3f}')
reflect = lambda c: CubicBezier(-c.start.conjugate(), -c.control1.conjugate(), -c.control2.conjugate(), -c.end.conjugate())

class Cluster:
    # tiles is a list of edge-lists, where each edge-list is a SVG path
    # keys is a 4-long list of (in_edge, out_edge) where in_edge.end = out_edge.start = key point
    # colors is a string with length equal to tiles
    # keys are reversed (both order and start/end) because first operation in step is to reverse them
    def __init__(self, edges, keys, color):
        self.tiles = [edges]
        self.keys = [(edges[k+1], edges[k]) for k in reversed(keys)]
        self.colors = color
    
    def step(self, color):
        # Reflect across y-axis (keeps layout of cluster clockwise)
        self.tiles = [[reflect(edge) for edge in tile] for tile in self.tiles]
        self.keys = [tuple(reflect(edge) for edge in reversed(pair)) for pair in reversed(self.keys)]

        seen = set(mark for tile in self.tiles for edge in tile for mark in marks(edge))
        tiles = []
        keys = []
        colors = ''

        # keys[0].in to keys[0].out (first new key point)
        I = bezier_isometry(self.keys[0][0], self.keys[0][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        keys += [pathkeys[1]]
        colors += ''.join(color for tile in add)

        # keys[1].in to pathkeys[3].out (no key point)
        I = bezier_isometry(self.keys[1][0], pathkeys[3][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        colors += ''.join(color for tile in add)

        # keys[0].in to pathkeys[2].out (no key point)
        I = bezier_isometry(self.keys[0][0], pathkeys[2][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        colors += ''.join(color for tile in add)

        # keys[1].in to pathkeys[3].out (second new key point)
        I = bezier_isometry(self.keys[1][0], pathkeys[3][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        keys += [pathkeys[2]]
        colors += ''.join(color for tile in add)

        # keys[1].in to pathkeys[3].out (no key point)
        I = bezier_isometry(self.keys[1][0], pathkeys[3][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        colors += ''.join(color for tile in add)

        # keys[0].in to pathkeys[2].out (third new key point)
        I = bezier_isometry(self.keys[0][0], pathkeys[2][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        keys += [pathkeys[1]]
        colors += ''.join(color for tile in add)

        # keys[1].in to pathkeys[3].out (fourth new key point)
        I = bezier_isometry(self.keys[1][0], pathkeys[3][1])
        add = [[bezier_apply(I, e) for e in path] for path in self.tiles]
        add = [[e for e in path if marks(e)[0] not in seen] for path in add]
        add = [a for a in add if len(a) > 0]
        pathkeys = [tuple(bezier_apply(I, edge) for edge in pair) for pair in self.keys]
        seen.update(set(mark for tile in add for edge in tile for mark in marks(edge)))
        tiles += add
        keys += [pathkeys[2]]
        colors += ''.join(color for tile in add)

        self.tiles += tiles
        self.keys = keys
        self.colors += colors

    # key has two edges, find the point common to both
    @classmethod
    def key_point(cls, key):
        a,b = key
        s = set([f'{zround(a.start,3):.3f}', f'{zround(a.end,3):.3f}'])
        return b.start if f'{zround(b.start,3):.3f}' in s else b.end

    def svg(self, path):
        # Set SVG dimensions so each tile is approximately 50mm at widest point
        user_unit = 38.6;
        pair = [(Path(segment), k) for tile, k in zip(self.tiles, self.colors) for segment in tile]
        paths = [item[0] for item in pair]
        colors = ''.join(item[1] for item in pair)
        stroke_widths = [1/user_unit for item in pair]
        nodes = [Cluster.key_point(k) for k in self.keys]
        node_radii = [3/user_unit for k in self.keys]

        minx = min(min(path.start.real, path.end.real) for path in paths) - 1
        miny = min(min(path.start.imag, path.end.imag) for path in paths) - 1
        maxx = max(max(path.start.real, path.end.real) for path in paths) + 1
        maxy = max(max(path.start.imag, path.end.imag) for path in paths) + 1

        viewbox = (minx, miny, maxx - minx, maxy - miny)
        dimensions = tuple(user_unit*i for i in viewbox[-2:])

        wsvg(paths=paths, filename=path, nodes=nodes, colors=colors, 
             stroke_widths=stroke_widths, node_radii=node_radii,
             viewbox=viewbox, dimensions=dimensions)

parser = argparse.ArgumentParser(
                    prog='algorithm-a',
                    description='Generate a tiling using a chiral aperiodic monotile',
                    epilog='epilog')
parser.add_argument('-g', '--generations', help='number of genrations to run algorithm A', default=0, type=int)
parser.add_argument('path', help='output path for SVG file')
options = parser.parse_args()

# Create initial tile
points = reduce(lambda a, i: a + [a[-1] + i], moves[1:-1], [moves[0]]) + [0j]
keys = [0,2,4,8]
edges = [bezier(points[i], points[i+1], i) for i in range(14)]
c = Cluster(edges, keys, 'r')

# Apply the algorithm to expand the cluster
colors = 'cmbkluy'
for i in range(options.generations):
    c.step(colors[i % len(colors)])

# Write out cluster
c.svg(options.path)

