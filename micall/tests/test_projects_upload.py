import pytest

from micall.utils.projects_upload import (
    fetch_pipeline_by_version,
    fetch_region_by_name,
    find_missing_scoring_regions,
)


class DummyProjectConfig:
    def __init__(self):
        self.config = {
            'projects': {
                'DemoProject': {
                    'regions': [
                        {'coordinate_region': 'AA1', 'seed_region_names': ['AA1-seed']},
                        {'coordinate_region': 'NUC1', 'seed_region_names': ['NUC1-seed']},
                    ]
                }
            },
            'regions': {
                'AA1': {'is_nucleotide': False},
                'NUC1': {'is_nucleotide': True},
            },
        }


class DummySession:
    def __init__(self, regions):
        self._regions = regions

    def get_json(self, path):
        assert path == '/lab_miseq_regions'
        return self._regions


class SequencedSession:
    def __init__(self, responses_by_path):
        self._responses_by_path = responses_by_path
        self._calls_by_path = {path: 0 for path in responses_by_path}

    def get_json(self, path):
        responses = self._responses_by_path[path]
        call_index = self._calls_by_path[path]
        self._calls_by_path[path] += 1
        if call_index >= len(responses):
            return responses[-1]
        return responses[call_index]


def test_find_missing_scoring_regions_ignores_nucleotide_regions():
    project_config = DummyProjectConfig()
    scoring_config = {
        'projects': {
            'DemoProject': {
                'regions': []
            }
        }
    }

    missing = find_missing_scoring_regions(project_config, scoring_config)

    assert missing == [('DemoProject', 'AA1')]


def test_fetch_region_by_name_returns_match():
    session = DummySession([
        {'name': 'RT', 'id': 1},
        {'name': 'CA', 'id': 2},
    ])

    region = fetch_region_by_name(session, 'CA')

    assert region == {'name': 'CA', 'id': 2}


def test_fetch_region_by_name_raises_if_missing():
    session = DummySession([{'name': 'RT', 'id': 1}])

    with pytest.raises(RuntimeError, match="Region 'CA' was not returned by QAI"):
        fetch_region_by_name(session, 'CA', max_attempts=2, wait_seconds=0)


def test_fetch_region_by_name_retries_until_visible():
    session = SequencedSession({
        '/lab_miseq_regions': [
            [{'name': 'RT', 'id': 1}],
            [{'name': 'CA', 'id': 2}],
        ]
    })

    region = fetch_region_by_name(session, 'CA', max_attempts=3, wait_seconds=0)

    assert region == {'name': 'CA', 'id': 2}


def test_fetch_pipeline_by_version_returns_match():
    session = SequencedSession({
        '/lab_miseq_pipelines?version=7.18': [
            [],
            [{'id': 17, 'version': '7.18'}],
        ]
    })

    pipeline = fetch_pipeline_by_version(session, '7.18', max_attempts=3, wait_seconds=0)

    assert pipeline == {'id': 17, 'version': '7.18'}


def test_fetch_pipeline_by_version_raises_if_missing():
    session = SequencedSession({
        '/lab_miseq_pipelines?version=7.18': [[], []]
    })

    with pytest.raises(RuntimeError, match="Pipeline '7.18' was not returned by QAI"):
        fetch_pipeline_by_version(session, '7.18', max_attempts=2, wait_seconds=0)
