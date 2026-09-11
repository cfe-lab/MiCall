import os
import tarfile

from micall.utils.micall_kive import create_coverage_maps_tar

CONTENTS = {
    "conseq-cov.png": b"fake-png-bytes-1",
    "nuc-cov.png": b"fake-png-bytes-2",
    "amino-cov.png": b"fake-png-bytes-3",
}


def make_maps_dir(base, names_in_order, mtime, mode):
    maps_dir = base / "maps"
    maps_dir.mkdir(parents=True)
    for name in names_in_order:
        path = maps_dir / name
        path.write_bytes(CONTENTS[name])
        os.utime(path, (mtime, mtime))
        os.chmod(path, mode)
    return maps_dir


def test_coverage_maps_tar_is_deterministic(tmp_path):
    maps1 = make_maps_dir(tmp_path / "run1",
                          sorted(CONTENTS),
                          mtime=1000000000,
                          mode=0o644)
    maps2 = make_maps_dir(tmp_path / "run2",
                          sorted(CONTENTS, reverse=True),
                          mtime=1700000000,
                          mode=0o600)
    tar1 = tmp_path / "maps1.tar"
    tar2 = tmp_path / "maps2.tar"

    create_coverage_maps_tar(str(tar1), str(maps1))
    create_coverage_maps_tar(str(tar2), str(maps2))

    assert tar1.read_bytes() == tar2.read_bytes()


def test_coverage_maps_tar_contents(tmp_path):
    maps_dir = make_maps_dir(tmp_path / "run",
                             sorted(CONTENTS, reverse=True),
                             mtime=1700000000,
                             mode=0o600)
    tar_path = tmp_path / "maps.tar"

    create_coverage_maps_tar(str(tar_path), str(maps_dir))

    expected_names = [f"coverage_maps/{name}" for name in sorted(CONTENTS)]
    with tarfile.open(str(tar_path), mode="r") as tar:
        members = tar.getmembers()
        assert [m.name for m in members] == expected_names
        for member in members:
            assert member.isfile()
            assert member.mtime == 0
            assert member.uid == 0
            assert member.gid == 0
            assert member.uname == ""
            assert member.gname == ""
            assert member.mode & 0o777 == 0o644
            with tar.extractfile(member) as f:
                assert f is not None
                base = member.name.split("/")[-1]
                assert f.read() == CONTENTS[base]


def test_coverage_maps_tar_hardlinks(tmp_path):
    payload = b"shared-bytes"
    maps1 = tmp_path / "run1" / "maps"
    maps1.mkdir(parents=True)
    for name in ("a.png", "b.png"):
        (maps1 / name).write_bytes(payload)
        os.utime(maps1 / name, (1000000000, 1000000000))
        os.chmod(maps1 / name, 0o644)
    maps2 = tmp_path / "run2" / "maps"
    maps2.mkdir(parents=True)
    (maps2 / "a.png").write_bytes(payload)
    os.link(maps2 / "a.png", maps2 / "b.png")
    os.utime(maps2 / "a.png", (1700000000, 1700000000))
    os.chmod(maps2 / "a.png", 0o600)
    tar1 = tmp_path / "maps1.tar"
    tar2 = tmp_path / "maps2.tar"

    create_coverage_maps_tar(str(tar1), str(maps1))
    create_coverage_maps_tar(str(tar2), str(maps2))

    assert tar1.read_bytes() == tar2.read_bytes()
