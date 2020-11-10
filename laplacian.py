import math

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import numpy as np
import meshio

mesh = meshio.read('assets/domain.msh')

vertices = mesh.points

(inner_border_tag, _) = mesh.field_data['inner-border']
(outer_border_tag, _) = mesh.field_data['outer-border']
(domain_tag, _) = mesh.field_data['domain']

line_tags = mesh.cell_data_dict['gmsh:physical']['line']
triangle_tags = mesh.cell_data_dict['gmsh:physical']['triangle']

lines = mesh.cells_dict['line']
triangles = mesh.cells_dict['triangle']

# We do not want to include the sea vertices
inactive = np.full(len(vertices), False)
#inactive[ np.ndarray.flatten(lines[line_tags == border_tag]) ] = True
inactive[ np.ndarray.flatten(lines[line_tags == inner_border_tag]) ] = True
inactive[ np.ndarray.flatten(lines[line_tags == outer_border_tag]) ] = True

# To map the index in de 'vertices' array to an index
# into our matrix we need to subtract all inactive vertices
# that have come before the vertex we're currently 'looking' at.
subtract = np.cumsum(inactive)

map_index = np.arange(len(vertices)) - subtract

def interiormatrix(v1, v2, v3):
    tangent1 = v3 - v2
    tangent2 = v1 - v3
    tangent3 = v2 - v1
    normal = np.cross( (v1-v3), (v2-v3) )
    area = 0.5 * np.linalg.norm(normal)
    normal = normal / np.linalg.norm(normal) # Normalize
    grad1 = np.cross(normal, tangent1) / (2 *area)
    grad2 = np.cross(normal, tangent2) / (2 *area)
    grad3 = np.cross(normal, tangent3) / (2 *area)
      
    return area * np.array([
        [np.dot(grad1,grad1), np.dot(grad1,grad2), np.dot(grad1,grad3)],
        [np.dot(grad2,grad1), np.dot(grad2,grad2), np.dot(grad2,grad3)],
        [np.dot(grad3,grad1), np.dot(grad3,grad2), np.dot(grad3,grad3)]])

def elementvector(fun ,v1, v2, v3):
    el_size = np.linalg.norm(np.cross((v1-v3), (v2-v3))) / 2
    return el_size * np.array([fun(v1), fun(v2), fun(v3)]) / 3

N = np.sum(np.logical_not(inactive))
S = lil_matrix((N, N), dtype=np.float64)
F = np.zeros(N, dtype=np.float64)

print('Constructing matrix...')

# Let's assemble the matrix.
for el in triangles:
    v1, v2, v3 = vertices[el[0]], vertices[el[1]], vertices[el[2]]
    
    m = interiormatrix(v1, v2, v3)
    b = elementvector(lambda _: 5, v1, v2, v3)
    
    for p in range(0, 3):
        if inactive[el[p]]:
            continue
        
        for q in range(0, 3):
            if inactive[el[q]]:
                continue
            
            S[ map_index[el[p]], map_index[el[q]] ] += m[p, q]
         
        F[ map_index[el[p]] ] += b[p]

print('Solving matrix...')

u = spsolve(S.tocsr(), F)

print('Visualizing...')

z = np.zeros(len(vertices))
z[ np.logical_not(inactive) ] = u
z[ inactive ] = 0

x = vertices[:, 0]
y = vertices[:, 1]

triang = mtri.Triangulation(x, y, triangles)
plt.tricontourf(mtri.Triangulation(x, y, triangles), z, cmap='viridis', levels=50)
plt.colorbar()
plt.show()













