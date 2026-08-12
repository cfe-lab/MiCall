import argparse
import contextlib

import pytest

from micall.utils import projects_upload
from micall.utils.projects_upload import (
    create_pipeline,
    derive_order_by,
    fetch_pipeline_by_version,
    fetch_region_by_name,
    find_missing_scoring_regions,
    find_pipeline_by_version,
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
        self._calls_by_path = dict.fromkeys(responses_by_path, 0)

    def get_json(self, path):
        responses = self._responses_by_path[path]
        call_index = self._calls_by_path[path]
        self._calls_by_path[path] += 1
        if call_index >= len(responses):
            return responses[-1]
        return responses[call_index]


class RecordingSession:
    def __init__(self, pipelines):
        self._pipelines = pipelines
        self.posted = []

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False

    def login(self, qai_path, qai_user, password):
        pass

    def get_json(self, path):
        if path == '/lab_miseq_pipelines':
            return self._pipelines
        if path in ('/lab_miseq_seed_groups', '/lab_miseq_regions', '/lab_miseq_projects'):
            return []
        raise AssertionError('Unexpected get_json path: ' + path)

    def post_json(self, path, data):
        self.posted.append((path, data))
        return {'id': 1, **data}


class WrappedProjectsSession(RecordingSession):
    def get_json(self, path):
        if path == '/lab_miseq_projects':
            return {'content': [{'name': 'DemoProject', 'id': 7}]}
        if path == '/lab_miseq_regions':
            return [{'name': 'NUC1', 'is_nucleotide': True, 'reference': '',
                     'seed_group_id': None, 'id': 5}]
        return super().get_json(path)


class NucleotideOnlyProjectConfig:
    def __init__(self):
        self.config = {
            'projects': {
                'DemoProject': {
                    'regions': [
                        {'coordinate_region': 'NUC1',
                         'seed_region_names': ['NUC1-seed']},
                    ]
                }
            },
            'regions': {
                'NUC1': {'is_nucleotide': True, 'reference': [], 'seed_group': None},
            },
        }

    @classmethod
    def loadDefault(cls):
        return cls()


class PostingSession:
    def __init__(self, empty_response=False):
        self.empty_response = empty_response
        self.posted = []

    def get_json(self, path):
        assert path.startswith('/lab_miseq_pipelines?version=')
        return [{'id': 2, 'version': path.split('version=', 1)[1]}]

    def post_json(self, path, data):
        self.posted.append((path, data))
        if self.empty_response:
            return None
        return {'id': 1, **data}


class EmptyProjectConfig:
    def __init__(self):
        self.config = {'projects': {}, 'regions': {}}

    @classmethod
    def loadDefault(cls):
        return cls()


class ScoringFileStub:
    def __init__(self, path):
        self._path = path

    @contextlib.contextmanager
    def path(self):
        yield self._path


def make_args(pipeline_version='7.18', order_by=None):
    return argparse.Namespace(
        quiet=False,
        verbose=False,
        debug=False,
        qai_server='http://test-server',
        qai_user='bob',
        qai_password='testing',
        pipeline_version=pipeline_version,
        order_by=order_by,
        update_sequences=False)


def stub_main(monkeypatch, tmp_path, session, args,
              project_config_class=EmptyProjectConfig):
    monkeypatch.setattr(projects_upload, 'parse_args', lambda: args)
    monkeypatch.setattr(projects_upload, 'ProjectConfig', project_config_class)
    scoring_path = tmp_path / 'project_scoring.json'
    scoring_path.write_text('{"projects": {}}')
    monkeypatch.setattr(projects_upload, 'ProjectsScoringFile',
                        lambda: ScoringFileStub(scoring_path))
    monkeypatch.setattr(projects_upload.qai_helper, 'Session', lambda: session)


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


def test_find_pipeline_by_version():
    pipelines = [{'version': '1.0', 'order_by': 100},
                 {'version': '2.0', 'order_by': 200}]

    assert find_pipeline_by_version(pipelines, '2.0') == pipelines[1]
    assert find_pipeline_by_version(pipelines, '3.0') is None


def test_derive_order_by_uses_max_plus_100():
    pipelines = [{'version': '1.0', 'order_by': 100},
                 {'version': '2.0', 'order_by': 200}]

    assert derive_order_by(pipelines) == 300


def test_derive_order_by_defaults_to_100():
    assert derive_order_by([]) == 100
    assert derive_order_by([{'version': '1.0', 'order_by': None}]) == 100
    assert derive_order_by([{'version': '1.0'}]) == 100


def test_derive_order_by_uses_override():
    pipelines = [{'version': '1.0', 'order_by': 100}]

    assert derive_order_by(pipelines, override=50) == 50


def test_create_pipeline_posts_order_by():
    session = PostingSession()

    pipeline = create_pipeline(session, '7.18', 300)

    assert pipeline == {'id': 1, 'version': '7.18', 'order_by': 300}
    assert session.posted == [('/lab_miseq_pipelines',
                               {'version': '7.18', 'order_by': 300})]


def test_create_pipeline_fetches_by_version_on_empty_response():
    session = PostingSession(empty_response=True)

    pipeline = create_pipeline(session, '7.18', 300)

    assert pipeline == {'id': 2, 'version': '7.18'}
    assert session.posted == [('/lab_miseq_pipelines',
                               {'version': '7.18', 'order_by': 300})]


def test_main_posts_derived_order_by(monkeypatch, tmp_path):
    session = RecordingSession([{'version': '1.0', 'order_by': 100},
                                {'version': '2.0', 'order_by': 200}])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18'))

    assert projects_upload.main() == 0

    assert session.posted[0] == ('/lab_miseq_pipelines',
                                 {'version': '7.18', 'order_by': 300})


def test_main_posts_order_by_100_when_no_pipelines(monkeypatch, tmp_path):
    session = RecordingSession([])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18'))

    assert projects_upload.main() == 0

    assert session.posted[0][1]['order_by'] == 100


def test_main_posts_order_by_100_when_no_order_by_values(monkeypatch, tmp_path):
    session = RecordingSession([{'version': '1.0'},
                                {'version': '2.0', 'order_by': None}])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18'))

    assert projects_upload.main() == 0

    assert session.posted[0][1]['order_by'] == 100


def test_main_uses_order_by_override(monkeypatch, tmp_path):
    session = RecordingSession([{'version': '1.0', 'order_by': 100}])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18', order_by=50))

    assert projects_upload.main() == 0

    assert session.posted[0][1]['order_by'] == 50


def test_main_pipeline_already_exists_returns_1(monkeypatch, tmp_path):
    session = RecordingSession([{'version': '7.18', 'order_by': 300}])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18'))

    assert projects_upload.main() == 1

    assert session.posted == []


def test_main_unwraps_wrapped_projects_response(monkeypatch, tmp_path):
    session = WrappedProjectsSession([])
    stub_main(monkeypatch, tmp_path, session, make_args('7.18'),
              project_config_class=NucleotideOnlyProjectConfig)

    assert projects_upload.main() == 0

    assert session.posted[0] == ('/lab_miseq_pipelines',
                                 {'version': '7.18', 'order_by': 100})
    assert session.posted[1] == ('/lab_miseq_project_versions',
                                 {'pipeline_id': 1, 'project_id': 7})
    assert not any(path == '/lab_miseq_projects' for path, _ in session.posted)
