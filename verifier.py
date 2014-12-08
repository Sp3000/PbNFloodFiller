import PIL

filepath = input("Input filepath: ")
im = Image.open(filepath)
width, height = im.size

remaining_pixels = {(i, j) for i in range(height) for j in range(width)}
flood_fills = []
colors = set()

def neighbours(pixel):
    results = []

    for neighbour in [(pixel[0]+1, pixel[1]), (pixel[0]-1, pixel[1]),
                      (pixel[0], pixel[1]+1), (pixel[0], pixel[1]-1)]:

        if 0 <= neighbour[0] < width and 0 <= neighbour[1] < height:
            results.append(neighbour)

    return results

while remaining_pixels:
    start_pixel = remaining_pixels.pop()
    col = im.getpixel(start_pixel)

    to_search = {start_pixel}
    searched = set()
    ff = {start_pixel}

    while to_search:
        pixel = to_search.pop()

        for n in neighbours(pixel):
            if n not in searched and n not in ff:
                

    
