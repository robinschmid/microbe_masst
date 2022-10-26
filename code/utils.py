from pathlib import Path
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def prepare_paths(file=None, files=None):
    if files is not None:
        for f in files:
            try:
                Path(f).parent.mkdir(parents=True, exist_ok=True)
            except Exception as e:
                logger.exception(e)
    if file is not None:
        try:
            Path(file).parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            logger.exception(e)
