from micall.utils.spring_beads import Wire, Bead


def test_no_target():
    wire = Wire()
    wire.append(Bead(start=0, end=10))
    expected_displays = [(0, 10)]

    wire.align()
    assert wire.displays == expected_displays


def test_one_target():
    wire = Wire()
    wire.append(Bead(start=0, end=10, target_start=100, target_end=110))
    expected_start_displays = [(0, 10)]
    expected_aligned_displays = [(100, 110)]

    start_displays = wire.displays
    wire.align()
    aligned_displays = wire.displays

    assert start_displays == expected_start_displays
    assert aligned_displays == expected_aligned_displays


def test_target_before_loose():
    wire = Wire()
    wire.append(Bead(start=0, end=10, target_start=100, target_end=110))
    wire.append(Bead(start=10, end=20))
    expected_displays = [(100, 110), (110, 120)]

    wire.align()
    displays = wire.displays

    assert displays == expected_displays


def test_target_after_loose():
    wire = Wire()
    wire.append(Bead(start=0, end=10))
    wire.append(Bead(start=10, end=20, target_start=110, target_end=120))
    expected_displays = [(100, 110), (110, 120)]

    wire.align()
    displays = wire.displays

    assert displays == expected_displays


def test_loose_between_targets():
    wire = Wire()
    wire.append(Bead(start=0, end=10, target_start=0, target_end=10))
    wire.append(Bead(start=10, end=20))
    wire.append(Bead(start=20, end=30, target_start=100, target_end=110))
    expected_displays = [(0, 10), (50, 60), (100, 110)]

    wire.align()
    displays = wire.displays

    assert displays == expected_displays


def test_swapped_targets():
    wire = Wire()
    wire.append(Bead(start=0, end=10, target_start=90, target_end=100))
    wire.append(Bead(start=10, end=20, target_start=0, target_end=10))
    expected_displays = [(40, 50), (50, 60)]

    wire.align()
    displays = wire.displays

    assert displays == expected_displays


def test_iter():
    wire = Wire()
    wire.append(Bead(start=0, end=10, target_start=0, target_end=10))
    wire.append(Bead(start=10, end=20, target_start=100, target_end=120))
    expected_display_starts = [0, 105]

    wire.align()
    display_starts = [bead.display_start for bead in wire]

    assert display_starts == expected_display_starts


def test_init():
    wire = Wire([Bead(start=0, end=10, target_start=0, target_end=10),
                 Bead(start=10, end=20, target_start=100, target_end=120)])
    expected_displays = [(0, 10), (105, 115)]

    wire.align()

    displays = wire.displays
    assert displays == expected_displays


def test_repr():
    wire = Wire((Bead(start=0, end=10, target_start=0, target_end=10),
                 Bead(start=10, end=20, target_start=100, target_end=120)))
    expected_repr = 'Wire([Bead(0, 10), Bead(10, 20)])'

    wire.align()
    r = repr(wire)

    assert r == expected_repr
