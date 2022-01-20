"""
Error handling wrapper around naive_warp.py
"""

from absl import app
from absl import flags
import subprocess

FLAGS = flags.FLAGS

flags.DEFINE_string("out_file", "-", "Output, - means stdout")
flags.DEFINE_string(
    "err_file", "-", "Where to write output for failures, - means stderr"
)
flags.DEFINE_bool("border", True, "Whether to add a border to the flag")


def main(argv):
    assert len(argv) == 2
    cmd = ["python", "../naive_warp.py", "--out_file", FLAGS.out_file]
    if not FLAGS.border:
        cmd.extend(["--border_size", 0])
    cmd.append(argv[1])
    cmd_result = subprocess.run(cmd, capture_output=True, text=True)
    if cmd_result.returncode != 0:
        result = [
            "",
            f"ERROR {argv[1]}",
            f"returncode {cmd_result.returncode}",
        ]
        if cmd_result.stdout.strip():
            result.extend(
                [
                    "stdout",
                    "======",
                    cmd_result.stdout.strip(),
                    "",
                ]
            )

        if cmd_result.stderr.strip():
            result.extend(
                [
                    "stderr",
                    "======",
                    cmd_result.stderr.strip(),
                    "",
                ]
            )

        result = "\n".join(result)
        if FLAGS.err_file == "-":
            print(result)
        else:
            with open(FLAGS.err_file, "a") as f:
                f.write(result)


if __name__ == "__main__":
    app.run(main)
