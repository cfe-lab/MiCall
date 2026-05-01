import unittest

from micall.monitor.qai_helper import Session


class DummyResponse:
    def __init__(self, text, json_value=None, status_code=200):
        self.text = text
        self._json_value = json_value
        self.status_code = status_code

    def raise_for_status(self):
        return None

    def json(self):
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
