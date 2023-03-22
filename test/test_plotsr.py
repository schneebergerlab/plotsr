def test_cli():
    from subprocess import Popen, PIPE
    p = Popen('plotsr -h'.split(), stdout=PIPE, stderr=PIPE)
    out = p.communicate()
    assert 'error' not in out[1].decode()
