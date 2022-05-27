import pytest
from lorenz import enkf


def test_enkf_1():
    enkf.test_enkf_harmonique()