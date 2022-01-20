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

flags.DEFINE_bool("fail_on_error", True, "If True, fail if any file fails to save")

# internal flags, typically client wouldn't change
flags.DEFINE_string("build_dir", "build/", "Where build runs.")
flags.DEFINE_bool("exec_ninja", True, "Whether to run ninja.")
flags.DEFINE_string("noto_dir", None, "Dir to resolve Noto paths against")
flags.mark_flag_as_required("noto_dir")


def build_dir() -> Path:
    return Path(FLAGS.build_dir).resolve()


def noto_dir() -> Path:
    return Path(FLAGS.noto_dir).resolve()


def preflight_out_dir() -> Path:
    return build_dir() / "preflight"


def waved_out_dir() -> Path:
    return build_dir() / "waved"


def rel(from_path: Path, to_path: Path) -> Path:
    # relative_to(A,B) doesn't like it if B doesn't start with A
    return Path(os.path.relpath(str(to_path.resolve()), str(from_path.resolve())))


def rel_build(path: Path) -> Path:
    return rel(build_dir(), path)


def main(_) -> None:
    build_dir().mkdir(exist_ok=True)
    preflight_out_dir().mkdir(exist_ok=True)
    waved_out_dir().mkdir(exist_ok=True)
    build_file = build_dir() / "build.ninja"
    err_file = (build_dir() / "warp_errors.log").resolve()
    with open(build_file, "w") as f:
        nw = ninja_syntax.Writer(f)

        if FLAGS.fail_on_error:
            util_path = rel_build(Path("naive_warp.py"))
            nw.rule(f"waveflag", f"python {util_path} --out_file $out $in")
            nw.newline()
            nw.rule(
                f"waveflag-no-border",
                f"python {util_path} --border_size 0 --out_file $out $in",
            )
            nw.newline()
        else:
            util_path = rel_build(Path("warp_runner.py"))
            nw.rule(
                f"waveflag",
                f"python {util_path} --err_file {err_file} --out_file $out $in",
            )
            nw.newline()
            nw.rule(
                f"waveflag-no-border",
                f"python {util_path} --noborder --err_file {err_file} --out_file $out $in",
            )
            nw.newline()

        util_path = rel_build(Path("flag_preflight.py"))
        nw.rule(f"preflight", f"python {util_path} --out_file $out $in")
        nw.newline()

        with open("wave_list.txt") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                src_svg, wave_name, *options = (p.strip() for p in line.split(" "))
                waveflag_rule = "waveflag"
                if options and not bool(int(options[0])):
                    waveflag_rule = "waveflag-no-border"

                nw.build(
                    str(preflight_out_dir() / src_svg),
                    "preflight",
                    str(noto_dir() / src_svg),
                )

                nw.build(
                    str(waved_out_dir() / wave_name),
                    waveflag_rule,
                    str(preflight_out_dir() / src_svg),
                )

    ninja_cmd = ["ninja", "-C", os.path.dirname(build_file)]
    if FLAGS.exec_ninja:
        print(" ".join(ninja_cmd))
        if err_file.is_file():
            err_file.unlink()
        subprocess.run(ninja_cmd, check=True)


if __name__ == "__main__":

    app.run(main)
