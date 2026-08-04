from collections.abc import Iterator, Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import (
    Literal,
    Protocol,
    TypeAlias,
    TypedDict,
)

ALLOWED_GROUPS: Sequence[str] = ['Everyone']
# noinspection PyArgumentList
PipelineType = Enum(
    'PipelineType',
    ['FILTER_QUALITY', 'MAIN', 'MIDI', 'RESISTANCE', 'MIXED_HCV_MAIN', 'MIXED_HCV_MIDI',
     'DENOVO_MAIN', 'DENOVO_MIDI', 'DENOVO_RESISTANCE', 'PROVIRAL'])

PIPELINE_GROUPS: Mapping[PipelineType, PipelineType] = {
    PipelineType.FILTER_QUALITY: PipelineType.FILTER_QUALITY,
    PipelineType.MAIN: PipelineType.MAIN,
    PipelineType.MIDI: PipelineType.MAIN,
    PipelineType.RESISTANCE: PipelineType.MAIN,
    PipelineType.MIXED_HCV_MAIN: PipelineType.MIXED_HCV_MAIN,
    PipelineType.MIXED_HCV_MIDI: PipelineType.MIXED_HCV_MAIN,
    PipelineType.DENOVO_MAIN: PipelineType.DENOVO_MAIN,
    PipelineType.DENOVO_MIDI: PipelineType.DENOVO_MAIN,
    PipelineType.DENOVO_RESISTANCE: PipelineType.DENOVO_MAIN,
    PipelineType.PROVIRAL: PipelineType.PROVIRAL
}

PIPELINE_GROUP_DEPENDENCIES: Mapping[PipelineType, PipelineType] = {
    PipelineType.PROVIRAL: PipelineType.DENOVO_MAIN
}

class ConfigInterface(Protocol):
    micall_main_pipeline_id: int | None
    mixed_hcv_pipeline_id: int | None
    denovo_main_pipeline_id: int | None
    micall_filter_quality_pipeline_id: int | None
    micall_resistance_pipeline_id: int | None
    proviral_pipeline_id: int | None
    max_active: int
    pipeline_version: str
    kive_user: str
    kive_password: str
    kive_server: str
    raw_data: Path


class Item(TypedDict, total=False):
    """General type for items in a Kive batch."""
    name: str
    groups_allowed: Sequence[str]


class RunDataset(TypedDict):
    """Type for individual dataset entries in a run's dataset list."""
    id: str
    url: str
    argument_name: str
    argument_type: str
    dataset: str
    dataset_purged: bool
    name: str
    groups_allowed: Sequence[str]


class RunCreationDataset(TypedDict):
    """Type for dataset entries when creating a new run."""
    argument: str
    dataset: str


State: TypeAlias = Literal["C", "F", "X", "R", "N", "L"]


class Run(TypedDict):
    """TypedDict for a Kive run object.

    Based on usage in kive_watcher.py and sample_watcher.py:
    - 'id': string identifier used as key in active_runs
    - 'state': run state ('C' for complete, 'F' for failed, etc.)
    - 'datasets': list of run dataset entries (added when run completes)
    """
    id: str
    state: State
    datasets: list[RunDataset]


class Batch(Protocol):
    @property
    def name(self) -> str: ...
    def __getitem__(self, key: str) -> RunDataset: ...
    def get(self, key: str, default: RunDataset | None = None) -> RunDataset | None: ...
    def __iter__(self) -> Iterator[str]:  ...
