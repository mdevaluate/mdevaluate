import logging

logger = logging.getLogger('mdevaluate')
logger.setLevel(logging.INFO)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
logger.addHandler(stream_handler)

formatter = logging.Formatter('{levelname:8}[{asctime}]:{funcName}: {message}', style='{')
stream_handler.setFormatter(formatter)


def setlevel(level, file=None):
    """
    Change the level of logging. If `file` is specified, logs are written to this file.
    """
    if isinstance(level, str):
        level = getattr(logging, level)
    logger.setLevel(level)
    if file is not None:
        handler = logging.FileHandler(file)
        handler.setLevel(level)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        stream_handler.setLevel(level)
