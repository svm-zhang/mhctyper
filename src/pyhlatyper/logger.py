import logging
from typing import cast


class Logger(logging.Logger):
    def __init__(self, name: str):
        super().__init__(name)
        self._configured = False

    def initialize(self, debug: bool = False) -> None:
        if self._configured:
            return

        self.setLevel(logging.DEBUG if debug else logging.INFO)

        # Create console handler with specific format
        if not self.hasHandlers():
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter(
                "[%(module)s.py] - [%(funcName)s:%(lineno)d] - [%(levelname)s]: %(message)s"
            )
            console_handler.setFormatter(formatter)
            self.addHandler(console_handler)
        self._configured = True


logging.setLoggerClass(Logger)
# logging.getLogger returns Logger type
# Although setLoggerClass(Logger) sets the logger object as Logger class,
# the getLogger method, according to typeshed, still returns Logger type
# Thats why LSP (e.g. pyright) complains about initialize is not an attribute
# https://github.com/microsoft/pyright/discussions/5418
logger = cast(Logger, logging.getLogger(__name__))
