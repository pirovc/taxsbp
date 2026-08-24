import os
import shlex
import subprocess
import sys

base_dir = os.path.dirname(__file__)
sample_input = f"{base_dir}/data/sample.tsv"
sample_tax = f"{base_dir}/data/sample.tax"

from taxsbp.taxsbp import main as taxsbp_main


def run(cmd):
    errcode = 1
    stdout = None
    stderr = None

    process = subprocess.Popen(
        shlex.split(cmd),
        shell=False,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    stdout, stderr = process.communicate()
    errcode = process.returncode

    return errcode, stdout, stderr


def test_cli_help():
    errcode, stdout, stderr = run("taxsbp --help")
    assert errcode == 0
    assert stdout.startswith("usage: taxsbp [-h]")
    assert stderr == ""


def test_cli_stdout():
    errcode, stdout, stderr = run(f"taxsbp -i {sample_input} -t {sample_tax} -s")
    assert errcode == 0
    assert len(stdout.rstrip().split("\n")) == 13
    assert len(stdout.rstrip().split("\t")) == 40
    assert len(stderr.rstrip().split(",")) == 17


def test_cli_file():
    output_file = f"{base_dir}/test_cli_file.tsv"
    errcode, stdout, stderr = run(
        f"taxsbp -i {sample_input} -t {sample_tax} -o {output_file}"
    )
    assert errcode == 0
    assert stdout == ""
    assert stderr == ""
    assert os.path.isfile(output_file)
    with open(output_file, "r") as of:
        assert len(of.readlines()) == 13
    os.remove(output_file)


def test_cli_error():
    errcode, _, _ = run("taxsbp -x")  # non existing argument
    assert errcode != 0


def test_main_help(capsys):
    sys.argv = ["taxsbp --help"]
    try:
        taxsbp_main()
    except SystemExit as e:
        exit_code = e.code if e.code is not None else 0
    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("usage: taxsbp [-h]")
    assert captured.err == ""


def test_main_stdout(capsys):
    sys.argv = ["taxsbp", "-i", sample_input, "-t", sample_tax, "--stats"]
    try:
        taxsbp_main()
    except SystemExit as e:
        exit_code = e.code if e.code is not None else 0
    assert exit_code == 0
    captured = capsys.readouterr()
    assert len(captured.out.rstrip().split("\n")) == 13
    assert len(captured.err.rstrip().split(",")) == 17


def test_main_file(capsys):
    output_file = f"{base_dir}/test_main_file.tsv"
    sys.argv = [
        "taxsbp",
        "-i",
        sample_input,
        "-t",
        sample_tax,
        "-o",
        output_file,
    ]
    try:
        taxsbp_main()
    except SystemExit as e:
        exit_code = e.code if e.code is not None else 0
    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
    with open(output_file, "r") as of:
        assert len(of.readlines()) == 13
    os.remove(output_file)


def test_main_error(capsys):
    exit_code = None
    sys.argv = ["taxsbp", "-x"]  # non existing argument
    try:
        taxsbp_main()
    except SystemExit as e:
        exit_code = e.code if e.code is not None else 0
    assert exit_code != 0
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("usage: taxsbp [-h]")
