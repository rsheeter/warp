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
SWITCH_TAG = f"{{{svgns()}}}switch"
FOREIGN_OBJECT_TAG = f"{{{svgns()}}}foreignObject"

FLAGS = flags.FLAGS

flags.DEFINE_string("out_file", "-", "Output, - means stdout")


# TODO: export as public function from picosvg?
def _swap_elements(swaps):
    for old_el, new_els in swaps:
        for new_el in reversed(new_els):
            old_el.addnext(new_el)
        parent = old_el.getparent()
        if parent is None:
            raise ValueError("Lost parent!")
        parent.remove(old_el)


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
    swaps = []
    for el in tree.getiterator("*"):

        # The AS.svg flag (American Samoa) uses unsupported <switch> element,
        # in turn containing an unsupported <foreignObject> element with some
        # Adobe Illustrator-specific metadata. This seem to be ignored by browsers.
        # Simply stripping the switch while keeping all its children (except
        # for the foreignObject one) fixes this specific flag. I am not sure if this
        # can be applied more generally to any switch in the wild, so for know
        # it's better to keep this hack in here, instead of in picosvg proper.
        if el.tag == SWITCH_TAG:
            swaps.append((el, [e for e in el if e.tag != FOREIGN_OBJECT_TAG]))
            continue

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

    _swap_elements(swaps)

    write_xml(FLAGS.out_file, tree, pretty=False)


if __name__ == "__main__":
    app.run(main)
