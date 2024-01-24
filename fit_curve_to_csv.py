"""
Fits a curve to a set of points taken from csv(s) of x,y pairs and outputs an svg of the result
"""

import csv
from pathlib import Path
import sys


from gg_fit_curve import fit_cubics


def read_csv(csv_file):
    with open(csv_file) as f:
        reader = csv.reader(f)
        rows = tuple(tuple(float(v) for v in row) for row in reader)
    return rows


def main():
    value_series = [read_csv(csv_file) for csv_file in sys.argv[1:]]

    max_sqr_err = 0.00001

    x_values = [x for series in value_series for (x,y) in series]
    y_values = [y for series in value_series for (x,y) in series]
    min_x = min(x_values)
    max_x = max(x_values)
    min_y = min(y_values)
    max_y = max(y_values)

    print(f"<svg viewBox=\"{min_x - 10} {min_y - 10}  {max_x - min_x + 20} {max_y - min_y + 20}\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">")
    for series in value_series:
        cubics = fit_cubics(series, max_sqr_err)
        print(f"fit_cubics max_sqr_err {max_sqr_err} cubics {len(cubics)}", file=sys.stderr)
        print(f"  <path fill=\"none\" stroke=\"black\" stroke-width=\"0.5\" d=\"")        
        for (i, (start, c0, c1, end)) in enumerate(cubics):
            if i == 0:
                print(f"          M {start[0]},{start[1]}")
            print(f"          C {c0[0]},{c0[1]} {c1[0]},{c1[1]} {end[0]},{end[1]}")
        print(f"\" />")

        for (x, y) in series:
            print(f"  <circle cx=\"{x}\" cy=\"{y}\" r=\"0.5\" fill=\"darkblue\" opacity=\"0.5\" />")
    print("</svg>")


if __name__ == "__main__":
    main()