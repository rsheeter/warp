"""
Noto region flags have some sharp edges; try to fix them before feeding to 
nanoemoji. Intended use is to run and commit the result to the noto-emoji repo.

Usage:
  python naive_warp.py noto-emoji/third_party/region-flags/svg/GB-WLS.svg --out_file GB-WLS-waved.svg
"""
from absl import app
from absl import flags

import enum
import math
from helpers import *
from lxml import etree
from picosvg.svg_meta import ntos, splitns, svgns, xlinkns, parse_view_box
import re

DEFAULT_PRECISION = 200
SHOULD_RENAME = {
    "store-width": "stroke-width",
}


class Dimension(enum.IntEnum):
    HORIZONTAL = 0
    VERTICAL = 1
    DIAGONAL = 2


# Currently a couple of flags (BR.svg and TW.svg) use percentages with <rect>,
# so we only deal with that for now. TODO remove once picosvg fully supports it
PERCENT_ATTRS = {
    f"{{{svgns()}}}rect": [
        (Dimension.HORIZONTAL, ("x", "width")),
        (Dimension.VERTICAL, ("y", "height")),
    ],
}

FLAGS = flags.FLAGS

flags.DEFINE_string("out_file", "-", "Output, - means stdout")


def _fix_pt(name, attrib):
    while match := re.search(
        r"(?:((?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+)))\s*pt", attrib[name]
    ):
        value = attrib[name]

        new_sz = ntos(round(float(match.group(1)) * 1.25, 2))
        new_value = value[: match.start()] + new_sz + value[match.end() :]
        attrib[name] = new_value


def _fix_percents(el, dimensions):
    if el.tag in PERCENT_ATTRS:
        for dim, attrs in PERCENT_ATTRS[el.tag]:
            scale = dimensions[dim.value]
            for attr in attrs:
                if attr in el.attrib and "%" in el.attrib[attr]:
                    n = float(el.attrib[attr].split("%")[0]) / 100 * scale
                    el.attrib[attr] = str(n)


def _parse_svg_dimensions(svg):
    if "viewBox" in svg.attrib:
        _, _, width, height = parse_view_box(svg.attrib["viewBox"])
    else:
        width = float(svg.attrib["width"])
        height = float(svg.attrib["height"])
    # cf. "normalized diagonal" formula at https://www.w3.org/TR/SVG11/coords.html#Units
    normalized_diagonal = math.hypot(width, height) / math.sqrt(2)
    return (width, height, normalized_diagonal)


def main(argv):
    tree = load_svg(argv, load_to=LoadTo.ETREE)

    dimensions = _parse_svg_dimensions(tree.getroot())

    el_to_rm = []
    for el in tree.getiterator("*"):
        # replace relative percentages with absolute numbers
        _fix_percents(el, dimensions)

        keys = list(el.attrib.keys())
        for name in keys:
            # find and eliminate pt's
            _fix_pt(name, el.attrib)

            # rename things
            if name in SHOULD_RENAME:
                new_name = SHOULD_RENAME[name]
                assert new_name not in keys
                el.attrib[new_name] = el.attrib[name]
                del el.attrib[name]

    write_xml(FLAGS.out_file, tree, pretty=False)


if __name__ == "__main__":
    app.run(main)
