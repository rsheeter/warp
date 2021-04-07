"""
Noto region flags have some sharp edges; try to fix them before feeding to 
nanoemoji. Intended use is to run and commit the result to the noto-emoji repo.

Usage:
  python naive_warp.py noto-emoji/third_party/region-flags/svg/GB-WLS.svg --out_file GB-WLS-waved.svg
"""
from absl import app
from absl import flags

import enum
from helpers import *
from lxml import etree
from picosvg.svg_meta import ntos, splitns, svgns, xlinkns
import re

DEFAULT_PRECISION = 200
SHOULD_RENAME = {
    "store-width": "stroke-width",
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


def main(argv):
    tree = load_svg(argv, load_to=LoadTo.ETREE)

    el_to_rm = []
    for el in tree.getiterator("*"):
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
