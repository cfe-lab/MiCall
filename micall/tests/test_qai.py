import micall.monitor.qai_helper as qai
import pytest
import csv
import os
import json
from pathlib import Path


@pytest.fixture
def conn():
    with qai.Session() as session:
        try:
            session.login('localhost:4567', 'qai', os.environ.get('QAI_PWD'))
        except Exception:
            return {
                'session': session,
                'cwd': Path(os.path.realpath(__file__)).parent
            }


def test_get(conn):
    table = conn['cwd'] / 'data' / 'proviral_sample.csv'
    with open(table, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in map(json.dumps, reader):
            conn.session.post_json('/proviral/create', row)
