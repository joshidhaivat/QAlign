import unittest

class ConvertTest(unittest.TestCase):

  def test_dummy(self):
    """
    usage: python -m unittest qalign.convert_test.ConvertTest.test_dummy
    """
    print('test_dummy')
    return

  def test_get_kmermap(self):
    from qalign.convert import get_kmermap

    path = None
    kmermap = get_kmermap(path)
    print(kmermap)

    return

  def test_get_qlevels(self):
    from qalign.convert import get_kmermap, get_qlevels 

    kmermap = get_kmermap()

    num_level_list = [2, 3, 4]
    for num_level in num_level_list: 
      qlevels = get_qlevels(kmermap, num_level)
      print('num_level=%d'%num_level)
      print(qlevels)

    return

  def test_convert_reads(self):
    import os
    from qalign.convert import convert 

    pa_dir = os.path.dirname(os.path.abspath(__file__))
    read_path = os.path.join(pa_dir, 'data', 'dummy_reads.fasta')
    test_dir = os.path.join(pa_dir, 'test')

    convert(read_path, test_dir, 2, rc=True, kmermap_path=None)

    return

if __name__ == '__main__':
  unittest.main()