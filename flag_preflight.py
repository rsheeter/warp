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
TEXT_TAG = f"{{{svgns()}}}text"
SWITCH_TAG = f"{{{svgns()}}}switch"
FOREIGN_OBJECT_TAG = f"{{{svgns()}}}foreignObject"


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
    svg_root = tree.getroot()

    has_text = svg_root.find(TEXT_TAG) is not None

    dimensions = _parse_svg_dimensions(svg_root)

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

            if name == "style":
                # picosvg doesn't support the 'display' attribute, and there's one flag
                # AC.svg (Ascension Island) which uses a style="display:inline" but that
                # is redundant since 'inline' is already the default value for that that
                # property and thus can be safely omitted.
                if "display:inline" in el.attrib[name]:
                    el.attrib[name] = el.attrib[name].replace("display:inline", "")

                # There's one flag TA.svg which has style="font-size:12x" everywhere
                # but doesn't even contain <text> elements to use that! Picosvg does
                # not currently support embedded text and complains about the unknown
                # font-size attribute, so we strip it below.
                if not has_text and "font-size" in el.attrib[name]:
                    el.attrib[name] = re.sub("font-size:[^;]+;?", "", el.attrib[name])

    _swap_elements(swaps)

    write_xml(FLAGS.out_file, tree, pretty=False)


if __name__ == "__main__":
    app.run(main)
