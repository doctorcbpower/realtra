"""Parameter files: TOML (native) and legacy Fortran namelists."""
import dataclasses
import sys
from pathlib import Path

import pytest

from galaxy_components import params as P
from galaxy_components.params import GROUP_FIELDS, Params, read_params, to_toml

DATA = Path(__file__).parent / "data"
EXAMPLES = Path(__file__).parent.parent / "examples"


def write(tmp_path, text, name="p.toml"):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_every_field_is_in_exactly_one_group():
    names = [n for g in GROUP_FIELDS.values() for n in g]
    assert sorted(names) == sorted(f.name for f in dataclasses.fields(Params))


@pytest.mark.parametrize("name", ["galaxy_params.txt", "nbody_ref_params.txt"])
def test_namelist_and_converted_toml_agree(tmp_path, name):
    nml = read_params(DATA / name)
    toml = write(tmp_path, to_toml(nml))
    assert read_params(toml) == nml
    assert read_params(write(tmp_path, to_toml(nml, all_fields=True), "all.toml")) == nml


def test_toml_defaults_lists_and_types(tmp_path):
    p = read_params(write(tmp_path, """
[halo]
m200 = 100            # an int is fine for a float
[nbody]
ncomp = 2
npart_comp = [1000, 500]
form_comp = ["NFW", "Hernquist"]
istaper = true
"""))
    assert p.m200 == 100.0 and isinstance(p.m200, float)
    assert p.npart_comp == [1000, 500, 0, 0, 0]                  # padded with defaults
    assert p.form_comp == ["nfw", "hernquist", "", "", ""]       # lower-cased, padded
    assert p.istaper is True and p.c200 == Params().c200          # untouched fields keep defaults


@pytest.mark.parametrize("text,fragment", [
    ("[halos]\nm200 = 1\n", "unknown table [halos]"),
    ("[halo]\nm_200 = 1\n", "unknown field 'm_200'"),
    ("[halo]\nndisc = 5\n", "belongs in [gas_disc]"),               # right name, wrong table
    ("[halo]\nm200 = \"big\"\n", "must be float"),
    ("[halo]\nm200 = true\n", "must be float"),                      # a bool is not a number
    ("[gas_disc]\nndisc = 1.5\n", "must be int"),
    ("[nbody]\nnpart_comp = 5\n", "must be a list"),
    ("[nbody]\nnpart_comp = [1, 2, 3, 4, 5, 6]\n", "at most 5"),
    ("[nbody]\nnpart_comp = [1.5]\n", "must be integers"),
    ("[nbody]\nform_comp = [1]\n", "must be strings"),
    ("[model]\ngalic = 1\n", "must be bool"),
    ("[halo\nm200 = 1\n", "invalid TOML"),
])
def test_toml_is_validated_strictly(tmp_path, text, fragment):
    with pytest.raises(ValueError, match=fragment.replace("[", r"\[").replace("]", r"\]")):
        read_params(write(tmp_path, text))


def test_suffix_selects_format_and_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_params(tmp_path / "nope.toml")
    # a .toml file is never parsed as a namelist, and a non-.toml file never as TOML
    nml_text = (DATA / "galaxy_params.txt").read_text()
    with pytest.raises(ValueError, match="invalid TOML"):
        read_params(write(tmp_path, nml_text, "looks_like_nml.toml"))
    assert read_params(write(tmp_path, nml_text, "x.nml")).m200 == 100.0


def test_namelist_without_f90nml_gives_a_helpful_error(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "f90nml", None)             # makes `import f90nml` fail
    with pytest.raises(ImportError, match="f90nml"):
        read_params(DATA / "galaxy_params.txt")
    # TOML needs nothing extra
    assert read_params(write(tmp_path, "[halo]\nm200 = 3.0\n")).m200 == 3.0


def test_example_toml_is_valid_and_describes_the_same_galaxy():
    toml = read_params(EXAMPLES / "galaxy_params.toml")           # strict: any typo would raise
    nml = read_params(EXAMPLES / "galaxy_params.txt")
    for name in ("m200", "c200", "mdisc", "fdisc", "mstar", "temp", "ndisc", "fbaryon", "nhalo",
                 "rtrunc", "form", "ncomp", "istaper", "halo_type", "galic"):
        assert getattr(toml, name) == getattr(nml, name), name
    assert toml.npart_comp[0] == nml.npart_comp[0] and toml.form_comp[0] == nml.form_comp[0]


def test_to_toml_only_writes_non_defaults():
    text = to_toml(Params(m200=50.0, ndisc=10))
    assert "m200 = 50.0" in text and "ndisc = 10" in text
    assert "c200" not in text and "[bulge]" not in text
    assert P.tomllib.loads(text) == {"halo": {"m200": 50.0}, "gas_disc": {"ndisc": 10}}
