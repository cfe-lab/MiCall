import json
import logging
import sys
import textwrap

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
        self.qai_user = qai_user

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
        full_path = self.qai_path + path
        response = method(
            full_path,
            data=json_data,
            headers=headers)

        if response.status_code >= 400:
            logger.info("Error response from QAI server %s", full_path)
            logger.info("Response status: %s", response.status_code)
            logger.info("Response body: %s", response.text)

        response.raise_for_status()
        if not response.text.strip():
            return None
        try:
            return response.json()
        except ValueError:
            redirects = ' -> '.join(
                f'{redirect.status_code} {redirect.headers.get("Location", "")}'
                for redirect in response.history)
            user = getattr(self, 'qai_user', None)
            user_hint = f' for user {user!r}' if user else ''
            print(f'QAI response body for {full_path}:', file=sys.stderr)
            print(textwrap.indent(response.text, '    '), file=sys.stderr)
            raise RuntimeError(
                f"QAI returned a non-JSON response{user_hint} for {full_path}: "
                f"{response.status_code} {response.reason}."
                + (f" Redirected: {redirects}." if redirects else "")
                + " The response body was printed above.") from None

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
