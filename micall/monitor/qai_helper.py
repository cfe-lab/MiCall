import json
import logging
import requests

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

    def _execute_json(self, method, path, data=None):
        """Execute a JSON request without retry logic.
        
        Retries should be handled by the calling code using retry_operation.
        """
        json_data = data and json.dumps(data)
        headers = {'Accept': 'application/json'}
        if json_data:
            headers['Content-Type'] = 'application/json'
            
        response = method(
            self.qai_path + path,
            data=json_data,
            headers=headers)
        response.raise_for_status()
        return response.json()

    def post_json(self, path, data):
        """ Post a JSON object to the web server, and return a JSON object.

        @param path the relative path to add to the qai_path used in login()
        @param data a JSON object that will be converted to a JSON string
        @return the response body, parsed as a JSON object
        """
        return self._execute_json(self.post, path, data)

    def get_json(self, path):
        """ Get a JSON object from the web server.

        @param path the relative path to add to QAI server path
        @return the response body, parsed as a JSON object
        """
        return self._execute_json(self.get, path)
