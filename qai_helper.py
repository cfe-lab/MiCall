import json
import logging
from random import Random
import requests
import time

import miseq_logging

logger = miseq_logging.init_logging_console_only(logging.INFO)

class Session(requests.Session):
    def login(self, qai_path, qai_user, password):
        """ Login to QAI before calling post_json or get_json.
        
        @raise RuntimeError: when the QAI server rejects the user and password.
        """
        self.qai_path = qai_path

        response = self.post(qai_path + "/account/login",
                             data={'user_login': qai_user,
                                   'user_password': password})
        if response.status_code == requests.codes.forbidden:  # @UndefinedVariable
            raise RuntimeError("Login failed for QAI user '{}'.".format(qai_user))

    def _retry_json(self, method, path, data=None, retries=3):
        json_data = data and json.dumps(data)
        retries_remaining = retries
        average_delay = 20
        while True:
            try:
                response = method(
                    self.qai_path + path,
                    data=json_data,
                    headers={'Content-Type': 'application/json',
                             'Accept': 'application/json'})
                response.raise_for_status()
                return response.json()
            except StandardError:
                if retries_remaining <= 0:
                    raise

                # ten minutes with some noise
                sleep_seconds = average_delay + Random().uniform(-10, 10)
                logger.warn(
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
        
        @param session an open HTTP session
        @param path the relative path to add to settings.qai_path
        @param retries: the number of times to retry the request before failing.
        @return the response body, parsed as a JSON object
        """
        return self._retry_json(self.get, path, retries=retries)
