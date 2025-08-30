import json
import logging
from random import Random
import requests
import time

logger = logging.getLogger('qai_helper')


class Session(requests.Session):
    def __init__(self) -> None:
        super().__init__()
        self.qai_path: str | None = None

    def login(self, qai_path: str, qai_user: str, password: str) -> None:
        """ Login to QAI before calling post_json or get_json.

        @raise RuntimeError: when the QAI server rejects the user and password.
        """
        self.qai_path = qai_path

        response = self.post(qai_path + "/account/login",
                             data={'user_login': qai_user,
                                   'user_password': password})
        response.raise_for_status()

    def _retry_json(self, method, path, data=None, retries=3):
        json_data = data and json.dumps(data)
        headers = {'Accept': 'application/json'}
        if json_data:
            headers['Content-Type'] = 'application/json'
        retries_remaining = retries
        average_delay = 20
        while True:
            # noinspection PyBroadException
            try:
                response = method(
                    self.qai_path + path,
                    data=json_data,
                    headers=headers)
                response.raise_for_status()
                return response.json()
            except Exception:
                if retries_remaining <= 0:
                    logger.error('JSON request failed for %s',
                                 path,
                                 exc_info=True)
                    raise

                # ten minutes with some noise
                sleep_seconds = average_delay + Random().uniform(-10, 10)
                logger.warning(
                    'JSON request failed. Sleeping for %ss before retry.',
                    sleep_seconds,
                    exc_info=True)
                time.sleep(sleep_seconds)
                retries_remaining -= 1
                average_delay += 600

    def post_json(self, path, data, retries=3):
        """ Post a JSON object to the web server, and return a JSON object.

        @param path the relative path to add to the qai_path used in login()
        @param data a JSON object that will be converted to a JSON string
        @param retries: the number of times to retry the request before failing.
        @return the response body, parsed as a JSON object
        """
        return self._retry_json(self.post, path, data, retries)

    def get_json(self, path, retries=3):
        """ Get a JSON object from the web server.

        @param path the relative path to add to QAI server path
        @param retries: the number of times to retry the request before failing.
        @return the response body, parsed as a JSON object
        """
        return self._retry_json(self.get, path, retries=retries)
