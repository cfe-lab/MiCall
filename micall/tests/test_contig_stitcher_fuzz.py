import pytest
import json
import os
from micall.utils.referencefull_contig_stitcher import (
    GenotypedContig,
    AlignedContig,
    stitch_consensus,
    stitch_contigs,
    drop_completely_covered,
    ReferencefullStitcherContext,
)
import micall.utils.referencefull_contig_stitcher as stitcher
import micall.utils.registry as registry
from micall.core.plot_contigs import build_stitcher_figure
from aligntools import CigarHit, Cigar, CigarActions
from typing import Dict, List
from collections import defaultdict


@pytest.fixture
def no_aligner(monkeypatch):
    monkeypatch.setattr("micall.utils.referencefull_contig_stitcher.align_to_reference", lambda x: [x])


@pytest.fixture(autouse=True)
def stitcher_context():
    stitcher.ReferencefullStitcherContext.set(stitcher.ReferencefullStitcherContext())
    registry.set(registry.Registry())


def read_contigs(line):
    array = json.loads(line)
    contig_descriptions = [obj["fields"] for obj in array if obj["type"] == "contig"]
    for description in contig_descriptions:
        start = description["start"]
        end = description["end"]
        name = description["name"]
        length = end - start + 1
        assert length > 0

        ref_seq = "A" * 1000  # it does not matter
        seq = "C" * 10 + "A" * length + "T" * 10
        query = GenotypedContig(
            name=name,
            seq=seq,
            ref_name="commonref",
            group_ref="commongroup",
            ref_seq=ref_seq,
            match_fraction=2 / 3,
            reads_count=None,
        )
        alignment = CigarHit(
            Cigar([(length, CigarActions.MATCH)]),
            q_st=20,
            q_ei=20 + length - 1,
            r_st=start,
            r_ei=end,
        )
        contig = AlignedContig.make(query=query, alignment=alignment, strand="forward")
        aidee = f"{start:03d}-{end:03d}"
        yield {"contig": contig, "id": aidee}


def get_case_descriptions():
    pwd = os.path.dirname(__file__)
    jsonfile = os.path.join(pwd, "data", "contig_stitcher_fuzz_nogaps.json")
    with open(jsonfile, "r", encoding="utf8") as reader:
        for line in reader:
            read = list(read_contigs(line))
            contigs = [x["contig"] for x in read]
            ids = [x["id"] for x in read]
            aidee = ",".join(ids)
            yield {"contigs": contigs, "id": aidee}


all_case_descriptions = list(get_case_descriptions())
all_case_ids = [x["id"] for x in all_case_descriptions]


@pytest.mark.parametrize("description", all_case_descriptions, ids=all_case_ids)
def test_contig_number_prop(no_aligner, description):
    contigs = description["contigs"]
    stitched = list(stitch_consensus(contigs))
    assert len(stitched) <= len(contigs)


@pytest.mark.parametrize("description", all_case_descriptions, ids=all_case_ids)
def test_contig_number_prop2(no_aligner, description):
    contigs = description["contigs"]
    consensus = list(stitch_consensus(contigs))
    stitched = list(stitch_contigs(contigs))
    uncovered = list(drop_completely_covered(contigs))
    assert len(consensus) <= len(stitched) <= len(uncovered) <= len(contigs)


def test_contig_number_prop2_existential():
    # This test is just to confirm that our cases cover all sub-actions.

    contig_sets = [x["contigs"] for x in all_case_descriptions]

    assert any(
        len(list(stitch_contigs(contigs))) > len(list(stitch_consensus(contigs)))
        for contigs in contig_sets
    )

    assert any(
        len(list(drop_completely_covered(contigs))) > len(list(stitch_contigs(contigs)))
        for contigs in contig_sets
    )

    assert any(
        len(list(contigs)) > len(list(drop_completely_covered(contigs)))
        for contigs in contig_sets
    )


def get_all_reference_positions(contigs: List[GenotypedContig]):
    ret: Dict[int, int] = defaultdict(lambda: 0)
    for contig in contigs:
        if isinstance(contig, AlignedContig):
            for i in contig.alignment.coordinate_mapping.ref_to_query.domain:
                ret[i] += 1

    return ret


@pytest.mark.parametrize("description", all_case_descriptions, ids=all_case_ids)
def test_stitching_intervals_prop(no_aligner, description):
    contigs = description["contigs"]
    stitched = list(stitch_contigs(contigs))
    initial_positions = get_all_reference_positions(contigs)
    stitched_positions = get_all_reference_positions(stitched)

    # Checks that no reference position has been lost, and no new positions "created"
    assert set(initial_positions.keys()) == set(stitched_positions.keys())

    # Checks that there are no overlaps between contigs
    assert all(v == 1 for (k, v) in stitched_positions.items())


@pytest.mark.parametrize("description", all_case_descriptions, ids=all_case_ids)
def test_visualizer_simple(no_aligner, description):
    contigs = description["contigs"]
    with ReferencefullStitcherContext.fresh() as ctx:
        list(stitch_consensus(contigs))
        assert len(ctx.events) >= len(contigs)
        figure = build_stitcher_figure(ctx.events)
        assert len(figure.elements) > len(contigs) + 1
