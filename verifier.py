from PIL import Image

filepath = "out.png"
im = Image.open(filepath)
width, height = im.size

ungrouped_pixels = {(i, j) for i in range(width) for j in range(height)}
groups = []
colors = set()

def neighbours(pixel):
    results = []

    for neighbour in [(pixel[0]+1, pixel[1]), (pixel[0]-1, pixel[1]),
                      (pixel[0], pixel[1]+1), (pixel[0], pixel[1]-1)]:

        if 0 <= neighbour[0] < width and 0 <= neighbour[1] < height:
            results.append(neighbour)

    return results

while ungrouped_pixels:
    start_pixel = ungrouped_pixels.pop()
    group_color = im.getpixel(start_pixel)

    to_search = [start_pixel]
    group = {start_pixel}

    while to_search:
        pixel = to_search.pop()

        for n in neighbours(pixel):
            if n in ungrouped_pixels and im.getpixel(n) == group_color:
                ungrouped_pixels.remove(n)
                group.add(n)
                to_search.append(n)

    groups.append(group)
    colors.add(group_color)

print "N = {}, P = {}".format(len(groups), len(colors))
