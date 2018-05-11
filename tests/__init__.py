import logging
import pathlib


tests_directory = pathlib.Path(__file__).parent
data_directory = tests_directory / 'data'
logging.getLogger().setLevel(51)
