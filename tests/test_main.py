# python
import pytest

from main import greet, main


def test_greet():
    assert greet("Alice") == "Hi, Alice"


def test_greet_empty():
    assert greet("") == "Hi, "


def test_main_default(capsys):
    rc = main([])
    captured = capsys.readouterr()
    assert rc == 0
    assert "Hi, PyCharm" in captured.out


def test_main_custom(capsys):
    rc = main(["-n", "Bob"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "Hi, Bob" in captured.out
