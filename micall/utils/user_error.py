
class UserError(RuntimeError):
    """
    Base class for all exceptions that are to be presented to a user.
    """

    def __init__(self, fmt: str, *fmt_args: object):
        self.fmt = fmt
        self.fmt_args = fmt_args
        self.code = 1
        super().__init__(fmt % fmt_args if fmt_args else fmt)
