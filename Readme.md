# Python3 Tiling Generator for a Chiral Aperiodic Monotile

![Cartoon image of laser cut tiles](https://github.com/vincentl/chiral-aperiodic-tile/blob/main/images/banner.jpeg)

This project contains a python3 implementation of the algorithm described 
in Appendix A of [arXiv:2305.17743 [math.CO]](https://doi.org/10.48550/arXiv.2305.17743)
and several sample output tilings. One application of the tiling files is
for laser cutting multiple tiles with minimal waste.

# What is an Chiral Aperiodic Monotile?

- **monotile** = single tile can tile the xy-plane
- **aperiodic** = the tiling pattern _never_ repeats
- **chiral** = there are incompatible left handed and right handed versions of the tile

Tilings (or tessellations) have been studied for a very long time, but this
generator is based on a discovery first published in 2023.  It is a solution
to the [einstein problem](https://en.wikipedia.org/wiki/Einstein_problem) (see Wikipedia).

# Usage

Algorithm A proceeds in generations. We call the starting tile generation 0. The number of 
tiles grows quickly, so generations much above 5 are expensive to compute and difficult
to display.

Pick a generation number (for example 2) and an output file (for example `demo.svg`). Then 
run

```bash
python3 algorithm-a.py --generations 2 demo.svg
```

# Laser Cutting Tiles

Since the tiles are aperiodic, there is not a simple way to arrange the tiles for
laser cutting. Simply spacing the tiles in a grid leaves significant waste between tiles.
The algorithm for producing a tiling does not produce a simple rectangular tiling. By
starting with a tiling with a large enough rectangular patch for the intended cut-stock,
use a program like [Inkscape](https://inkscape.org/) to remove the paths beyond the
desired cut area.

The base tile has 14 paths. Each generation packs seven copies of the previous supertile
to form a new supertile. The python3 code is careful to no duplicate the edge paths
where supertiles are glued together. The end result is that the tiling files are composed
on many short paths. This may present challenges when cutting the paths depending
on the laser, software, and settings.

## Laser Cutting Example

![200mm square plywood sheet half yellow and half blue just after laser cutting](https://github.com/vincentl/chiral-aperiodic-tile/blob/main/images/laser-bed.jpeg)

![individual laser cut tiles ](https://github.com/vincentl/chiral-aperiodic-tile/blob/main/images/tiles.jpeg)

![laser cut tiles reassembled into a tiling](https://github.com/vincentl/chiral-aperiodic-tile/blob/main/images/tiling.jpeg)

# Sample Tiles

The `examples` directory contains sample tiles. The colors are used to group paths based
on which generation produced the path - 0 is red, 1 is cyan, etc. The four red nodes are
the anchor nodes used in Algorithm A to orient the current generation supertile in
producing the next generation supertile.

# Optional Tile Label

Experimental code to add a text label to each tile is available through command line options. A 
font name, size, and the label string are all required to activate the feature.

```bash
python3 algorithm-a.py --generations 2 --label 'Example' --font Helvetica --font-size 5.1mm label.svg
```
# Optional Control of Tile Border Colors

By default, the python3 code uses a different border color for each generation. To change the list of 
colors or force a single color, the command line option `--color` takes a string of color characters and
cycles through the list of colors for each generation. By making the list one long, a single color is used.

```bash
python3 algorithm-a.py --generations 3 --color r red.svg
```


# License

Python3 Tiling Generator for a Chiral Aperiodic Monotile © 2024 by Vincent
Lucarelli is licensed under [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/).
