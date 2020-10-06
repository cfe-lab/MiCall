import micall.monitor.qai_helper as qai
import pytest
import csv
import os
import json
import logging
from pathlib import Path


@pytest.fixture
def conn():
    with qai.Session() as session:
        try:
            session.login('http://localhost:4567', 'qai',
                          os.environ.get('QAI_PWD'))
        except Exception as e:
            logging.error(
                'SESSION FAILED TO LOG IN, did you set the password env var? "export QAI_PWD="'
            )
            logging.error(e)
        return {
            'session': session,
            'cwd': Path(os.path.realpath(__file__)).parent
        }


def test_create(conn):
    table = conn['cwd'] / 'data' / 'proviral_sample.csv'
    with open(table, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            conn['session'].post_json('/proviral/create', row)
