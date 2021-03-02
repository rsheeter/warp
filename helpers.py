import enum
from lxml import etree
from picosvg.svg import SVG


class LoadTo(enum.Enum):
    PICOSVG = enum.auto()
    ETREE = enum.auto()


def load_svg(argv, load_to: LoadTo = LoadTo.PICOSVG):
    assert len(argv) <= 2
    try:
        input_file = argv[1]
    except IndexError:
        input_file = None

    if load_to == LoadTo.PICOSVG:
        ldr = SVG
    elif load_to == LoadTo.ETREE:
        ldr = etree
    else:
        raise ValueError(f"Invalid load_to {load_to}")

    if input_file:
        svg = ldr.parse(input_file)
    else:
        svg = ldr.fromstring(sys.stdin.read())
    return svg


def _reduce_text(text):
    text = text.strip() if text else None
    return text if text else None


def reduce_whitespace(tree):
    # lxml really likes to retain whitespace
    for e in tree.iter("*"):
        e.text = _reduce_text(e.text)
        e.tail = _reduce_text(e.tail)


def write_xml(out_file, tree, pretty=True):
    out_content = etree.tostring(tree, pretty_print=pretty).decode("utf-8")
    if out_file == "-":
        print(out_content)
    else:
        with open(out_file, "w") as f:
            f.write(out_content)