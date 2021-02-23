"""Writes a ninja build to zopflipng and emplace all emoji *except* flags.

TODO: this can be ported to be a proper google3 tool
      IFF android_fonts equivalent is available [emoji_metadata?]
"""

from ninja import ninja_syntax
import os
from pathlib import Path
import subprocess
from typing import Sequence

from absl import app
from absl import flags


FLAGS = flags.FLAGS

# internal flags, typically client wouldn't change
flags.DEFINE_string("build_dir", "build/", "Where build runs.")
flags.DEFINE_bool("exec_ninja", True, "Whether to run ninja.")
flags.DEFINE_string("noto_dir", None, "Dir to resolve Noto paths against")


def build_dir() -> Path:
    return Path(FLAGS.build_dir).resolve()


def noto_dir() -> Path:
    return Path(FLAGS.noto_dir).resolve()


def out_dir() -> Path:
    return Path("waved")


def rel(from_path: Path, to_path: Path) -> Path:
    # relative_to(A,B) doesn't like it if B doesn't start with A
    return Path(os.path.relpath(str(to_path.resolve()), str(from_path.resolve())))


def rel_build(path: Path) -> Path:
    return rel(build_dir(), path)


def main(_) -> None:
  build_dir().mkdir(exist_ok=True)
  out_dir().mkdir(exist_ok=True)
  build_file = build_dir() / "build.ninja"
  with open(build_file, "w") as f:
    nw = ninja_syntax.Writer(f)

    util_path = rel_build(Path("naive_warp.py"))
    nw.rule(
        f"waveflag",
        f"python {util_path} --curve_order 2 --out_file $out $in"
    )
    nw.newline()

    with open("wave_list.txt") as f:
      for line in f:
        src_svg, wave_name = line.split(" ")
        nw.build(str(out_dir() / wave_name.strip()), "waveflag", str(noto_dir() / src_svg.strip()))

  ninja_cmd = ["ninja", "-C", os.path.dirname(build_file)]
  if FLAGS.exec_ninja:
      print(" ".join(ninja_cmd))
      subprocess.run(ninja_cmd, check=True)


if __name__ == '__main__':
  app.run(main)
