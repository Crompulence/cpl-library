from cplpy import CPL, cart_create, run_test, prepare_config, parametrize_file, CPL_VAR_TYPES

TESTS_DIR_NAMES = ["initialisation", "mapping"]


def exec_tests(test="all"):
    import pytest
    import os
    test_path = os.path.dirname(os.path.realpath(__file__))
    test_path = os.path.join(test_path, "test")
    print test_path
    if test != "all":
        test_path = os.path.join(test_path, test)
    pytest.main(["-v", test_path])

def get_test_dir():
    import os
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "test")
