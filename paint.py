from __future__ import division
from PIL import Image
import random, math
from collections import Counter

TOLERANCE = 0.2
INFILE = "mona.jpg"
OUTFILE = "out.png"
P = 30
N = 500

# https://stackoverflow.com/questions/2353211/hsl-to-rgb-color-conversion
def hsl(rgb):
    r, g, b = rgb
    R, G, B = r/255, g/255, b/255
    M = max(R, G, B)
    m = min(R, G, B)
    h = s = l = (M+m)/2

    if M == m:
        h = s = 0

    else:
        d = M - m
        s = d/(2 - M - m) if l > 0.5 else d/(M + m)

        if M == R:
            h = (G - B)/d + (6 if G < B else 0)
        elif M == G:
            h = (B - R)/d + 2
        elif M == B:
            h = (R - G)/d + 4

        h /= 6

    return (h, s, l)

def rgb(hsl):
    h, s, l = hsl
    
    if s == 0:
        r = g = b = l

    else:
        def hue_to_rgb(p, q, t):
            t = t % 1
            if t < 1/6: return p + (q-p)*6*t
            if t < 1/2: return q
            if t < 2/3: return p + (q-p)*(2/3-t)*6
            return p

        q = l*(1+s) if l<0.5 else l+s-l*s
        p = 2*l - q

        r = hue_to_rgb(p, q, h+1/3)
        g = hue_to_rgb(p, q, h)
        b = hue_to_rgb(p, q, h-1/3)

    return (int(round(r*255)), int(round(g*255)), int(round(b*255)))

im = Image.open(INFILE)
hsl_pixels = {}
width, height = im.size
out_im = Image.new("RGB", im.size)

for i in xrange(width):
    for j in xrange(height):
        hsl_pixels[(i, j)] = hsl(im.getpixel((i, j)))

print "Stage 1: HSL conversion complete"

def d(hsl1, hsl2):
    return ((min(abs(hsl1[0]-hsl2[0]), abs(hsl2[0]-hsl1[0]))/180)**2 + abs(hsl1[1]-hsl2[1])**2 + abs(hsl1[2]-hsl2[2])**2)**0.5
    # return max(abs(hsl1[0]-hsl2[0])/360, abs(hsl1[1]-hsl2[1]), abs(hsl1[2]-hsl2[2]))

def neighbours(pix):
    x = []

    for i in [(pix[0]+1, pix[1]), (pix[0]-1, pix[1]), (pix[0], pix[1]+1), (pix[0], pix[1]-1)]:
        if 0 <= i[0] < width and 0 <= i[1] < height:
            x.append(i)

    return x

def flood_fill(xy):
    to_fill = {xy}
    added = set()
    searched = set()
    col = hsl_pixels[xy]

    while to_fill:
        pix = to_fill.pop()
        pix_col = hsl_pixels[pix]

        if d(col, pix_col) < TOLERANCE:
            added.add(pix)
            unplaced.remove(pix)

            for i in neighbours(pix):
                if i in unplaced and i not in added and i not in searched:
                    to_fill.add(i)

        else:
            searched.add(pix)

    return added

# https://en.wikipedia.org/wiki/Directional_statistics
def mean(hsl_cols):
    mean_s = sum(x[1] for x in hsl_cols)/len(hsl_cols)
    mean_l = sum(x[2] for x in hsl_cols)/len(hsl_cols)

    s = sum(map(lambda x: math.sin(x[0]), hsl_cols))/len(hsl_cols)
    c = sum(map(lambda x: math.cos(x[0]), hsl_cols))/len(hsl_cols)

    if s >= 0 and c >= 0:
        mean_h = math.atan(s/c)
    elif c < 0:
        mean_h = math.atan(s/c) + 180
    elif s < 0 and c >= 0:
        mean_h = math.atan(s/c) + 360

    return (mean_h, mean_s, mean_l)

def mean_from_set(s):
    return rgb(mean([hsl_pixels[x] for x in s]))

def setnum(p):
    for i, s in enumerate(sets):
        if p in s:
            return i
        
sets = []
unplaced = {(i, j) for i in xrange(width) for j in xrange(height)}

while unplaced:
    point = unplaced.pop()
    unplaced.add(point)
    ff = flood_fill(point)
    sets.append(ff)

print "Stage 2: Flood fill set initiation complete"

while len(sets) > N:
    if len(sets) % 1000 == 0:
        print "%d sets remaining..." % len(sets)
        
    min_index = None
    min_elems = None

    for i in range(len(sets)):
        if min_elems == None or len(sets[i]) < min_elems:
            min_index = i
            min_elems = len(sets[i])

    merge_set = sets.pop(min_index)
    boundary = set()

    for p in merge_set:
        for neighbour in neighbours(p):
            if neighbour not in merge_set:
                boundary.add(neighbour)

    most_adjacent = Counter()

    for p in boundary:
        most_adjacent[setnum(p)] += 1

    merge_to_index = most_adjacent.most_common()[0][0]
    merge_to = sets[merge_to_index]
    merge_to |= merge_set
    sets[merge_to_index] = merge_to

print "Stage 3: Set merging complete"

print sorted(hsl_pixels[x][0] for x in sets[0])

set_cols = {}

for i,s in enumerate(sets):
    set_cols[i] = hsl(mean_from_set(s))

while len(set(set_cols.values())) > P:
    candidates = []
    
    for i in range(len(sets)):
        for j in range(i+1, len(sets)):
            if set_cols[i] != set_cols[j]:
                candidates.append((d(set_cols[i], set_cols[j]), i, j))

    # Merge with biased mean based on set size
    _, a, b = min(candidates)
    col1, col2 = set_cols[a], set_cols[b]
    size1, size2 = len(sets[a]), len(sets[b])

    new_s = (col1[1]*size1 + col2[1]*size2)/(size1 + size2)
    new_l = (col1[2]*size1 + col2[2]*size2)/(size1 + size2)

    s = (math.sin(col1[0])*size1 + math.sin(col2[0])*size2)/(size1 + size2)
    c = (math.cos(col1[0])*size1 + math.cos(col2[0])*size2)/(size1 + size2)

    if s >= 0 and c >= 0:
        new_h = math.atan(s/c)
    elif c < 0:
        new_h = math.atan(s/c) + 180
    elif s < 0 and c >= 0:
        new_h = math.atan(s/c) + 360
        
    new_col = (new_h, new_s, new_l)

    for i in range(len(sets)):
        if set_cols[i] in [col1, col2]:
            set_cols[i] = new_col

print "Stage 4: Colour merging complete"

for i in range(len(sets)):
    for p in sets[i]:
        out_im.putpixel(p, rgb(set_cols[i]))

out_im.save(OUTFILE)
print "Done!"
    
