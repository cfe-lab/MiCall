import pytest

from micall.core.project_config import ProjectConfig
from micall.utils.primer_tracker import PrimerTracker


@pytest.fixture(scope='session', name='projects')
def load_projects() -> ProjectConfig:
    return ProjectConfig.loadDefault()


def test_hcv1a(projects):
    conseq = projects.getReference('HCV-1a')
    tracker = PrimerTracker(conseq, 'HCV-1a')

    is_ignored_before = tracker.is_ignored(8244)
    is_ignored_within = tracker.is_ignored(8245)
    is_ignored_after = tracker.is_ignored(8282)

    assert not is_ignored_before
    assert is_ignored_within
    assert not is_ignored_after


def test_hiv(projects):
    conseq = projects.getReference('HIV1-B-FR-K03455-seed')
    tracker = PrimerTracker(conseq, 'HIV1-B-FR-K03455-seed')

    is_ignored = tracker.is_ignored(8400)

    assert not is_ignored


def test_deletion(projects):
    ref = projects.getReference('HCV-1a')
    conseq = ref[:8260] + '---' + ref[8263:]
    tracker = PrimerTracker(conseq, 'HCV-1a')

    is_ignored_before = tracker.is_ignored(8244)
    is_ignored_within = tracker.is_ignored(8281)
    is_ignored_after = tracker.is_ignored(8282)

    assert not is_ignored_before
    assert is_ignored_within
    assert not is_ignored_after
