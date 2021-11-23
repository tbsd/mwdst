import sys
import os
from os.path import exists
from PIL import Image, ImageDraw
import time

def graph_to_img(path):
    lines = {}
    with open(path) as fp:
        lines = fp.readlines()
    size = int(lines[0].split()[-1])
    points = [(int(line.split()[0]), int(line.split()[1])) for line in lines[1:]]

    minX = minY = maxX = maxY = 0;
    for p in points:
        if (minX > p[0]):
            minX = p[0]
        if (minY > p[1]):
            minY = p[1]
        if (maxX < p[0]):
            maxX = p[0]
        if (maxY < p[1]):
            maxY = p[1]

    dX = -minX
    dY = -minY

    im = Image.new(mode="RGB", size = (maxX + dX + 1, maxY + dY + 1))

    for p in points:
        im.putpixel((p[0] + dX, p[1] + dY), (255, 255, 255))

    im.save("{}.png".format(path))

def draw_solution(path, edges_path):
    if not exists(edges_path):
        return
    lines = {}
    with open(path) as fp:
        lines = fp.readlines()
    size = int(lines[0].split()[-1])
    points = [(int(line.split()[0]), int(line.split()[1])) for line in lines[1:]]

    minX = minY = maxX = maxY = 0;
    for p in points:
        if (minX > p[0]):
            minX = p[0]
        if (minY > p[1]):
            minY = p[1]
        if (maxX < p[0]):
            maxX = p[0]
        if (maxY < p[1]):
            maxY = p[1]

    dX = -minX
    dY = -minY

    im = Image.new(mode="RGB", size = (maxX + dX + 1, maxY + dY + 1))

    for p in points:
        im.putpixel((p[0] + dX, p[1] + dY), (255, 255, 255))

    lines = {}
    with open(edges_path) as fp:
        lines = fp.readlines()
    edges = [(int(line.split()[1]) - 1, int(line.split()[2]) - 1) for line in lines[2:]]
    draw = ImageDraw.Draw(im)
    for e in edges:
        source = points[e[0]]
        target = points[e[1]]
        edge = [source, target]
        draw.line(edge, fill = "red", width = 0)

    im.save("{}_solution.png".format(path))

def draw_with_edges(path, edges_path):
    if not exists(edges_path):
        return
    lines = {}
    with open(path) as fp:
        lines = fp.readlines()
    size = int(lines[0].split()[-1])
    points = [(int(line.split()[0]) * 10, int(line.split()[1]) * 10) for line in lines[1:]]

    minX = minY = maxX = maxY = 0;
    for p in points:
        if (minX > p[0]):
            minX = p[0]
        if (minY > p[1]):
            minY = p[1]
        if (maxX < p[0]):
            maxX = p[0]
        if (maxY < p[1]):
            maxY = p[1]

    dX = -minX
    dY = -minY

    im = Image.new(mode="RGB", size = (maxX + dX + 1, maxY + dY + 1))

    
    edges = [(int(line.split()[0]) - 1, int(line.split()[1]) - 1) for line in lines[1:]]
    draw = ImageDraw.Draw(im)
    #  for e in edges:
        #  source = points[e[0]]
        #  target = points[e[1]]
        #  edge = [source, target]
        #  draw.line(edge, fill = "blue", width = 0)

    lines = {}
    with open(edges_path) as fp:
        lines = fp.readlines()
    edges = [(int(line.split()[1]) - 1, int(line.split()[2]) - 1) for line in lines[2:]]
    for e in edges:
        source = points[e[0]]
        target = points[e[1]]
        edge = [source, target]
        draw.line(edge, fill = "red", width = 0)


        #  source = points[16]
        #  target = points[10]
        #  edge = [source, target]
        #  draw.line(edge, fill = "blue", width = 0)
        #  source = points[42]
        #  target = points[16]
        #  edge = [source, target]
        #  draw.line(edge, fill = "green", width = 0)

    for p in points:
        im.putpixel((p[0] + dX, p[1] + dY), (255, 255, 255))


    im.save("{}_with_edges.png".format(edges_path))

def draw_with_edges2(path, edges_path):
    if not exists(edges_path):
        return
    lines = {}
    with open(path) as fp:
        lines = fp.readlines()
    size = int(lines[0].split()[-1])
    points = [(int(line.split()[0]), int(line.split()[1])) for line in lines[1:]]

    minX = minY = maxX = maxY = 0;
    for p in points:
        if (minX > p[0]):
            minX = p[0]
        if (minY > p[1]):
            minY = p[1]
        if (maxX < p[0]):
            maxX = p[0]
        if (maxY < p[1]):
            maxY = p[1]

    dX = -minX
    dY = -minY

    im = Image.new(mode="RGB", size = (maxX + dX + 1, maxY + dY + 1))

    
    draw = ImageDraw.Draw(im)
    for p1 in points:
        for p2 in points:
            edge = [p1, p2]
            draw.line(edge, fill = "blue", width = 0)

    lines = {}
    with open(edges_path) as fp:
        lines = fp.readlines()
    edges = [(int(line.split()[1]) - 1, int(line.split()[2]) - 1) for line in lines[2:]]
    for e in edges:
        source = points[e[0]]
        target = points[e[1]]
        edge = [source, target]
        draw.line(edge, fill = "red", width = 0)
        if (e[0] == 1 and e[1] == 2 or e[0] == 2 and e[1] == 1):
            draw.line(edge, fill = "blue", width = 0)
        if (e[0] == 1 and e[1] == 7) or e[0] == 7 and e[1] == 1:
            draw.line(edge, fill = "green", width = 0)

    for p in points:
        im.putpixel((p[0] + dX, p[1] + dY), (255, 255, 255))


    im.save("{}_with_edges.png".format(path))

#  directory = 'tmpPics'
#  for filename in os.listdir(directory):
    #  f = os.path.join(directory, filename)
    #  if os.path.isfile(f):
        #  draw_with_edges("Taxicab_512.txt", f)
#  draw_with_edges("Taxicab_512.txt", "Kurbatov_512.txt")
#  exit(0)

#  while True:

#  graph_to_img("Taxicab_64.txt")
#  graph_to_img("Taxicab_128.txt")
#  graph_to_img("Taxicab_512.txt")
#  graph_to_img("Taxicab_2048.txt")
#  graph_to_img("Taxicab_4096.txt")
#
#  draw_solution("Taxicab_64.txt", "Kurbatov_64.txt")
#  draw_solution("Taxicab_128.txt", "Kurbatov_128.txt")
#  draw_solution("Taxicab_512.txt", "Kurbatov_512.txt")
#  draw_solution("Taxicab_2048.txt", "Kurbatov_2048.txt")
#  draw_solution("Taxicab_4096.txt", "Kurbatov_4096.txt")
#
draw_with_edges("Taxicab_64.txt", "Kurbatov_64.txt")
draw_with_edges("Taxicab_128.txt", "Kurbatov_128.txt")
draw_with_edges("Taxicab_512.txt", "Kurbatov_512.txt")
draw_with_edges("Taxicab_2048.txt", "Kurbatov_2048.txt")
draw_with_edges("Taxicab_4096.txt", "Kurbatov_4096.txt")
    #  exit(0)

#  draw_with_edges2("Taxicab_64.txt", "Kurbatov_64.txt")
#  draw_with_edges2("Taxicab_128.txt", "Kurbatov_128.txt")
#  draw_with_edges2("Taxicab_512.txt", "Kurbatov_512.txt")
#  draw_with_edges2("Taxicab_2048.txt", "Kurbatov_2048.txt")
#  draw_with_edges2("Taxicab_4096.txt", "Kurbatov_4096.txt")
    #  time.sleep(1)
