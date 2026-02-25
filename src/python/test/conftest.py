import pytest


def pytest_addoption(parser):
    parser.addoption("--filepath", action="store", default="")


def pytest_configure(config):
    config.addinivalue_line("markers", "filepath: mark test as in need of file path")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--filepath"):
        # --filepath given in cli: do not skip tests
        return
    skip_filepath = pytest.mark.skip(reason="needs --filepath option to run")
    for item in items:
        if "filepath" in item.keywords:
            item.add_marker(skip_filepath)
