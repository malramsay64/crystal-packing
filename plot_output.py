#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import sys
import matplotlib.pyplot as plt 
import re

remove_chars = "[" + re.escape("{}[]()") + "]"

def main(filename):
    fig, ax = plt.subplots()
    with open(filename) as src:
        for line in src:
            if line.startswith("Atom"):
                line = re.sub(remove_chars, '', line[4:])
                x, y, r, colour =[element.strip() for element in line.strip().split(',')]
                circle = plt.Circle((float(x), float(y)), float(r), facecolor=colour,
                    edgecolor=colour, linewidth=3, alpha=0.8)
                ax.add_patch(circle)
            elif line.startswith("Line"):
                line = re.sub(remove_chars, '', line[4:])
                x0, y0, x1, y1, colour =[element.strip() for element in line.strip().split(',')]
                ax.plot([float(x0), float(x1)], [float(y0), float(y1)], c=colour)
            else:
                line = re.sub(remove_chars, '', line)
                x0, y0, x1, y1, colour =[element.strip() for element in line.strip().split(',')]
                ax.plot([float(x0), float(x1)], [float(y0), float(y1)], c=colour)
                
    ax.set_aspect(1)
    ax.plot(0, 0)

    fig.show()

    plt.savefig(filename.with_suffix('.pdf'))

if __name__ == "__main__":
    if len(sys.argv) == 2:
        try:
            filename = Path(argv[1])
        except ValueError:
            raise ValueError("Invalid filename")
    else:
        filename = Path("test.txt")

    if not filename.exists():
        raise ValueError(f"File '{filename}' does not exist")

    main(filename)

