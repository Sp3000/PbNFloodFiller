from __future__ import division
from PIL import Image
import random, math, time
from collections import Counter, defaultdict, namedtuple

"""
Configure settings here
"""

INFILE = "nebula.jpg"
OUTFILE = "out.png"
P = 30
N = 500

FLOOD_FILL_TOLERANCE = 10
CLOSE_cell_TOLERANCE = 5
SMALL_CELL_THRESHOLD = 10
K_MEANS_TRIALS = 20
FILTER_RADIUS = 5

"""
Color conversion functions
"""

# http://www.easyrgb.com/?X=MATH
def rgb2xyz(rgb):
    r, g, b = rgb
    r /= 255
    g /= 255
    b /= 255

    r = ((r + 0.055)/1.055)**2.4 if r > 0.04045 else r/12.92
    g = ((g + 0.055)/1.055)**2.4 if g > 0.04045 else g/12.92
    b = ((b + 0.055)/1.055)**2.4 if b > 0.04045 else b/12.92

    r *= 100
    g *= 100
    b *= 100

    x = r*0.4124 + g*0.3576 + b*0.1805
    y = r*0.2126 + g*0.7152 + b*0.0722
    z = r*0.0193 + g*0.1192 + b*0.9505

    return (x, y, z)

def xyz2lab(xyz):
    x, y, z = xyz
    x /= 95.047
    y /= 100
    z /= 108.883

    x = x**(1/3) if x > 0.008856 else 7.787*x + 16/116
    y = y**(1/3) if y > 0.008856 else 7.787*y + 16/116
    z = z**(1/3) if z > 0.008856 else 7.787*z + 16/116

    L = 116*y - 16
    a = 500*(x - y)
    b = 200*(y - z)

    return (L, a, b)

def rgb2lab(rgb):
    return xyz2lab(rgb2xyz(rgb))

def lab2xyz(lab):
    L, a, b = lab
    y = (L + 16)/116
    x = a/500 + y
    z = y - b/200

    y = y**3 if y**3 > 0.008856 else (y - 16/116)/7.787
    x = x**3 if x**3 > 0.008856 else (x - 16/116)/7.787
    z = z**3 if z**3 > 0.008856 else (z - 16/116)/7.787

    x *= 95.047
    y *= 100
    z *= 108.883

    return (x, y, z)

def xyz2rgb(xyz):
    x, y, z = xyz
    x /= 100
    y /= 100
    z /= 100

    r = x* 3.2406 + y*-1.5372 + z*-0.4986
    g = x*-0.9689 + y* 1.8758 + z* 0.0415
    b = x* 0.0557 + y*-0.2040 + z* 1.0570

    r = 1.055 * (r**(1/2.4)) - 0.055 if r > 0.0031308 else 12.92*r
    g = 1.055 * (g**(1/2.4)) - 0.055 if g > 0.0031308 else 12.92*g
    b = 1.055 * (b**(1/2.4)) - 0.055 if b > 0.0031308 else 12.92*b

    r *= 255
    g *= 255
    b *= 255

    return (r, g, b)

def lab2rgb(lab):
    rgb = xyz2rgb(lab2xyz(lab))
    return (int(round(rgb[0])), int(round(rgb[1])), int(round(rgb[2])))

"""
Stage 1: Read in image and convert to CIELAB
"""

im = Image.open(INFILE)
width, height = im.size

def make_pixlab_map(im):
    width, height = im.size
    pixlab_map = {}

    for i in xrange(width):
        for j in xrange(height):
            pixlab_map[(i, j)] = rgb2lab(im.getpixel((i, j)))

    return pixlab_map

pixlab_map = make_pixlab_map(im)

print "Stage 1: CIELAB conversion complete"

"""
Stage 2: Partitioning the image into like-colored cells using flood fill
"""

def d(color1, color2):
    return (abs(color1[0]-color2[0])**2 + abs(color1[1]-color2[1])**2 + abs(color1[2]-color2[2])**2)**.5

def neighbours(pixel):
    results = []
 
    for neighbour in [(pixel[0]+1, pixel[1]), (pixel[0]-1, pixel[1]),
                      (pixel[0], pixel[1]+1), (pixel[0], pixel[1]-1)]:
 
        if 0 <= neighbour[0] < width and 0 <= neighbour[1] < height:
            results.append(neighbour)
 
    return results

def flood_fill(start_pixel):
    to_search = {start_pixel}
    cell = set()
    searched = set()
    start_color = pixlab_map[start_pixel]

    while to_search:
        pixel = to_search.pop()

        if d(start_color, pixlab_map[pixel]) < FLOOD_FILL_TOLERANCE:
            cell.add(pixel)
            unplaced_pixels.remove(pixel)

            for n in neighbours(pixel):
                if n in unplaced_pixels and n not in cell and n not in searched:
                    to_search.add(n)

        else:
            searched.add(pixel)

    return cell

# These two maps are inverses, pixel/s <-> cell number
cell_sets = {}
pixcell_map = {}
unplaced_pixels = {(i, j) for i in xrange(width) for j in xrange(height)}

while unplaced_pixels:
    start_pixel = unplaced_pixels.pop()
    unplaced_pixels.add(start_pixel)
    cell = flood_fill(start_pixel)

    cellnum = len(cell_sets)
    cell_sets[cellnum] = cell

    for pixel in cell:
        pixcell_map[pixel] = cellnum

print "Stage 2: Flood fill partitioning complete, %d cells" % len(cell_sets)

"""
Stage 3: Merge cells with less than a specified threshold amount of pixels to reduce the number of cells
         Also good for getting rid of some noise
"""

for i in range(1, SMALL_CELL_THRESHOLD):
    small_cells = []

    for cellnum in cell_sets:
        if len(cell_sets[cellnum]) <= i:
            small_cells.append(cellnum)

    for cellnum in small_cells:
        neighbour_cells = []

        for cell in cell_sets[cellnum]:
            for n in neighbours(cell):
                neighbour_reg = pixcell_map[n]

                if neighbour_reg != cellnum:
                    neighbour_cells.append(neighbour_reg)

        closest_cell = max(neighbour_cells, key=neighbour_cells.count)

        for cell in cell_sets[cellnum]:
            pixcell_map[cell] = closest_cell

        if len(cell_sets[closest_cell]) <= i:
            small_cells.remove(closest_cell)

        cell_sets[closest_cell] |= cell_sets[cellnum]
        del cell_sets[cellnum]

        if len(cell_sets) <= N:
            break

    if len(cell_sets) <= N:
            break
        
print "Stage 3: Small cell merging complete, %d cells" % len(cell_sets)

"""
Stage 4: Close color merging
"""

def mean_color(cell):
    L_sum = 0
    a_sum = 0
    b_sum = 0

    for pixel in cell:
        L, a, b = pixlab_map[pixel]
        L_sum += L
        a_sum += a
        b_sum += b

    return L_sum/len(cell), a_sum/len(cell), b_sum/len(cell)

cell_means = {}

for cellnum in cell_sets:
    cell_means[cellnum] = mean_color(cell_sets[cellnum])

neighbour_graph = defaultdict(set)

for i in range(width):
    for j in range(height):
        pixel = (i, j)
        cell = pixcell_map[pixel]

        for n in neighbours(pixel):
            neighbour_cell = pixcell_map[n]

            if neighbour_cell != cell:
                neighbour_graph[cell].add(neighbour_cell)
                neighbour_graph[neighbour_cell].add(cell)

print "Stage 4a: Neighbour graph built"

def merge_cells(merge_from, merge_to):
    merge_from_cell = cell_sets[merge_from]

    for pixel in merge_from_cell:
        pixcell_map[pixel] = merge_to

    del cell_sets[merge_from]
    del cell_means[merge_from]

    neighbour_graph[merge_to] |= neighbour_graph[merge_from]
    neighbour_graph[merge_to].remove(merge_to)

    for n in neighbour_graph[merge_from]:
        neighbour_graph[n].remove(merge_from)

        if n != merge_to:
            neighbour_graph[n].add(merge_to)

    del neighbour_graph[merge_from]

    cell_sets[merge_to] |= merge_from_cell
    cell_means[merge_to] = mean_color(cell_sets[merge_to])

# Go through the cells from largest to smallest. Keep replenishing the list while we can still merge.
last_time = time.time()
to_search = sorted(cell_sets.keys(), key=lambda x:len(cell_sets[x]), reverse=True)
full_list = True

while len(cell_sets) > N and to_search:
    if time.time() - last_time > 15:
        last_time = time.time()
        print "Close color merging... (%d cells remaining)" % len(cell_sets)

    while to_search:
        cellnum = to_search.pop()
        close_cells = []

        for neighbour_cellnum in neighbour_graph[cellnum]:
            if d(cell_means[cellnum], cell_means[neighbour_cellnum]) < CLOSE_cell_TOLERANCE:
                close_cells.append(neighbour_cellnum)

        if close_cells:
            for neighbour_cellnum in close_cells:
                merge_cells(neighbour_cellnum, cellnum)

                if neighbour_cellnum in to_search:
                    to_search.remove(neighbour_cellnum)
    
            break

    if full_list == True:
        if to_search:
            full_list = False

    else:
        if not to_search:
            to_search = sorted(cell_sets.keys(), key=lambda x:len(cell_sets[x]), reverse=True)
            full_list = True

print "Stage 4b: Close color merging complete, %d cells" % len(cell_sets)

"""
Stage 5: N-merging - merge until <= N cells
         Want to merge either 1) small cells or 2) cells close in color
"""

# Weight score between neighbouring cells by 1) size of cell and 2) color difference
def score(cell1, cell2):
    return d(cell_means[cell1], cell_means[cell2]) * len(cell_sets[cell1])**.5

neighbour_scores = {}

for cellnum in cell_sets:
    for n in neighbour_graph[cellnum]:
        neighbour_scores[(n, cellnum)] = score(n, cellnum)

last_time = time.time()

while len(cell_sets) > N:
    if time.time() - last_time > 15:
        last_time = time.time()
        print "N-merging... (%d cells remaining)" % len(cell_sets)
        
    merge_from, merge_to = min(neighbour_scores, key=lambda x: neighbour_scores[x])

    for n in neighbour_graph[merge_from]:
        del neighbour_scores[(merge_from, n)]
        del neighbour_scores[(n, merge_from)]
        
    merge_cells(merge_from, merge_to)

    for n in neighbour_graph[merge_to]:
        neighbour_scores[(n, merge_to)] = score(n, merge_to)
        neighbour_scores[(merge_to, n)] = score(merge_to, n)
        
print "Stage 5: N-merging complete, %d cells" % len(cell_sets)

"""
Stage 6: P merging - use k-means
"""

def form_clusters(centroids):
    clusters = defaultdict(set)

    for cellnum in cell_sets:
        scores = []
        
        for centroid in centroids:
            scores.append((d(centroid, cell_means[cellnum]), centroid))

        scores.sort()
        clusters[scores[0][1]].add(cellnum)

    return clusters

def calculate_centroid(cluster):
    L_sum = 0
    a_sum = 0
    b_sum = 0
    
    for cellnum in cluster:
        color = cell_means[cellnum]

        L_sum += color[0]
        a_sum += color[1]
        b_sum += color[2]

    return (L_sum/len(cluster), a_sum/len(cluster), b_sum/len(cluster))

def db_index(clusters):
    # Davies-Bouldin index

    scatter = {}

    for centroid, cluster in clusters.items():
        scatter_score = 0

        for cellnum in cluster:
            scatter_score += d(cell_means[cellnum], centroid)

        scatter_score /= len(cluster)
        scatter[centroid] = scatter_score**2 # Mean squared distance

    index = 0

    for ci, cluster in clusters.items():
        dist_scores = []
        
        for cj in clusters:
            if ci != cj:
                dist_scores.append((scatter[ci] + scatter[cj])/d(ci, cj))

        index += max(dist_scores)

    return index

best_clusters = None
best_index = None

for i in range(K_MEANS_TRIALS):  
    centroids = {cell_means[cellnum] for cellnum in random.sample(cell_sets, P)}
    converged = False

    while not converged:
        clusters = form_clusters(centroids)
        new_centroids = {calculate_centroid(cluster) for cluster in clusters.values()}

        if centroids == new_centroids:
            converged = True

        centroids = new_centroids

    index = db_index(clusters)

    if best_index is None or index < best_index:
        best_index = index
        best_clusters = clusters

cell_colors = {}

for centroid, cluster in best_clusters.items():
    for cellnum in cluster:
        cell_colors[cellnum] = centroid

print "Stage 6: P-merging complete"

"""
Stage 7: Gaussian smoothing
"""

def G(dx, dy):
    return math.e**(-(dx**2 + dy**2)/FILTER_RADIUS)

FilterPoint = namedtuple("FilterPoint", ["dx", "dy", "g"])
g_filter = []

for dx in range(-FILTER_RADIUS, FILTER_RADIUS+1):
    for dy in range(-FILTER_RADIUS, FILTER_RADIUS+1):
        g_filter.append(FilterPoint(dx, dy, G(dx, dy))) 

def closest_color(c, pixel):
    possible_colors = set()

    for n in neighbours(pixel):
        possible_colors.add(cell_colors[pixcell_map[n]])
        
    return min(possible_colors, key=lambda x: d(x, c))

out_im = Image.new("RGB", im.size)

for i in range(width):
    for j in range(height):
        pixel = (i, j)
        L_sum = 0
        a_sum = 0
        b_sum = 0
        g_sum = 0

        for filter_point in g_filter:
            nearby_pixel = (i + filter_point.dx, j + filter_point.dy)

            if 0 <= nearby_pixel[0] < width and 0 <= nearby_pixel[1] < height:
                nearby_color = cell_colors[pixcell_map[nearby_pixel]]
                
                L_sum += filter_point.g * nearby_color[0]
                a_sum += filter_point.g * nearby_color[1]
                b_sum += filter_point.g * nearby_color[2]
                g_sum += filter_point.g

        avg_color = (L_sum/g_sum, a_sum/g_sum, b_sum/g_sum)
        out_im.putpixel(pixel, lab2rgb(closest_color(avg_color, pixel)))

out_im.save("test.png")
    

"""
Stage last: Output the image!
"""

out_im = Image.new("RGB", im.size)

for cellnum in cell_sets:
    color = lab2rgb(cell_colors[cellnum])

    for pixel in cell_sets[cellnum]:
        out_im.putpixel(pixel, color)

out_im.save(OUTFILE)

print "Done!"
