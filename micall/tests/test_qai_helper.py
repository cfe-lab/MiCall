import unittest

from micall.monitor.qai_helper import Session


class DummyResponse:
    def __init__(self, text, json_value=None, status_code=200, headers=None, history=None):
        self.text = text
        self._json_value = json_value
        self.status_code = status_code
        self.reason = 'OK'
        self.headers = headers or {}
        self.history = history or []

    def raise_for_status(self):
        return None

    def json(self):
        if self._json_value is None:
            raise ValueError('No JSON object could be decoded')
        return self._json_value


class QaiHelperSessionTest(unittest.TestCase):
    def test_execute_json_returns_none_for_empty_success(self):
        session = Session()
        session.qai_path = 'http://example.invalid'

        def fake_method(url, data, headers):
            return DummyResponse('   \n')

        response = session._execute_json(fake_method, '/lab_miseq_regions', {'name': 'CA'})

        self.assertIsNone(response)

    def test_execute_json_returns_parsed_payload_when_present(self):
        session = Session()
        session.qai_path = 'http://example.invalid'

        def fake_method(url, data, headers):
            return DummyResponse('{"id": 17}', {'id': 17})

        response = session._execute_json(fake_method, '/lab_miseq_regions', {'name': 'CA'})

        self.assertEqual({'id': 17}, response)

    def test_execute_json_raises_clear_error_on_non_json_response(self):
        session = Session()
        session.qai_path = 'http://example.invalid'
        session.qai_user = 'bob'
        history = [
            DummyResponse('', status_code=303, headers={'Location': 'http://example.invalid/'}),
            DummyResponse(
                '',
                status_code=302,
                headers={'Location': 'http://example.invalid/qcs_start/labtech'}),
        ]

        def fake_method(url, data, headers):
            return DummyResponse('<html>Start Page</html>', status_code=200, history=history)

        with self.assertRaisesRegex(
                RuntimeError,
                r"non-JSON response for user 'bob' for "
                r'http://example\.invalid/lab_miseq_pipelines: 200 OK\.'
                r' Redirected: 303 http://example\.invalid/ -> '
                r'302 http://example\.invalid/qcs_start/labtech\.'
                r' The QAI user may be missing the group required for this operation') as cm:
            session._execute_json(fake_method, '/lab_miseq_pipelines', {'version': '0-dev'})

        self.assertIn("user 'bob'", str(cm.exception))
        self.assertIn('/qcs_start/labtech', str(cm.exception))
