import unittest
import pathlib
import ismap_v2
import read_grouping
import mapping_to_ref
import mapping_to_query
import os
import create_output
from Bio import SeqIO

#TODO: have tests that use already created fastq files/bam files etc to speed up testing

class TestGetSeqs(unittest.TestCase):

    def setUp(self):
        self.query_single = ["/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta"]
        self.query_list = ["/Users/jane/Desktop/ismap_v2/queries/IS26.fasta",
                           "/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta"]
        self.query_multi = ["/Users/jane/Desktop/ismap_v2/queries/two_ISqueries.fasta"]

        self.genome_single = ["/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk"]
        self.genome_multi = ["/Users/jane/Desktop/ismap_v2/refs/two_genomes.gbk"]
        self.genome_list = ["/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk",
                            "/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk"]

    def test_get_seqs_01(self):
        # check that when we give valid files, we return a list object
        self.assertIsInstance(ismap_v2.get_sequences(self.query_single, 'fasta'), list)
        self.assertIsInstance(ismap_v2.get_sequences(self.query_list, 'fasta'), list)
        self.assertIsInstance(ismap_v2.get_sequences(self.query_multi, 'fasta'), list)

        self.assertIsInstance(ismap_v2.get_sequences(self.genome_single, 'genbank'), list)
        self.assertIsInstance(ismap_v2.get_sequences(self.genome_list, 'genbank'), list)
        self.assertIsInstance(ismap_v2.get_sequences(self.genome_multi, 'genbank'), list)

    def test_get_seqs_02(self):
        # check that if the list is empty, we print a useful error
        with self.assertRaises(ismap_v2.NoSeqError):
            ismap_v2.get_sequences([], 'fasta')


class TestCreate_typing_output(unittest.TestCase):


    def setUp(self):
        self.filenames = {'intersect':'/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP010781.1_intersect.bed',
                          'closest': '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP010781.1_closest.bed'}
        self.ref_gbk_obj = SeqIO.read('/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk', 'genbank')

    def test_create_typing_output(self):
        hit_list = create_output.create_typing_output(self.filenames, self.ref_gbk_obj)


'''
class TestMapToISQuery(unittest.TestCase):

    def setUp(self):
        # set up reads
        self.test_reads = [pathlib.Path('/Users/jane/Desktop/ismap_v2/reads/9262_1#29_1.fastq.gz'),
                      pathlib.Path('/Users/jane/Desktop/ismap_v2/reads/9262_1#29_2.fastq.gz')]
        self.read_groups = read_grouping.group_reads(self.test_reads)
        self.sample = self.read_groups.paired[0]

        # set up IS query
        self.query_single = ['/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta']
        self.query_records = ismap_v2.get_sequences(self.query_single, 'fasta')
        self.is_query = self.query_records[0]

        # set up output
        self.output = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29'

    def test_map_to_is_query_01(self):

        test_left, test_right= ismap_v2.map_to_is_query(self.sample, self.is_query, self.output)

        self.assertTrue(os.path.isfile(test_left))
        self.assertTrue(os.path.isfile(test_right))

        # verify right number

        # verify correct reads

        #self.assertEqual(ismap_v2.map_to_is_query(self.sample, self.is_query, self.output), ['/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/9262_1#29_ISAba1_left_final.fastq',
                                             #'/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/9262_1#29_ISAba1_right_final.fastq'])
'''
class TestSetOutputFilenames(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/tmp'
        self.prefix = '9262_1#29'
        self.query = 'ISAba1'
        self.out_dir = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/'

    def test_set_output_filenames(self):
        pass

class TestRefMapping(unittest.TestCase):

    def setUp(self):

         self.left_flanking = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_ISAba1_left_final.fastq'
         self.right_flanking = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_ISAba1_right_final.fastq'

         ref_seqs = ismap_v2.get_sequences(['/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk'], 'genbank')
         self.ref = ref_seqs[0]
         self.sample = '9262_1#29'
         self.tmp_folder = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/tmp'
         self.out_folder = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/'


    def test_map_to_ref_seq_01(self):

        test_filenames = mapping_to_ref.map_to_ref_seq(self.ref, self.sample, self.left_flanking,
                                                                            self.right_flanking, self.tmp_folder, self.out_folder, '1')
        test_left_sorted = test_filenames['left_merged_bed']
        test_right_sorted = test_filenames['right_merged_bed']

        self.assertTrue(os.path.isfile(test_left_sorted))
        self.assertTrue(os.path.isfile(test_right_sorted))

class TestCreateBedFiles(unittest.TestCase):

    def setUp(self):
        pass

        self.ref_name = 'CP010781.1'
        self.sample = '9262_1#29'
        self.tmp_folder = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/tmp'
        self.out_folder = '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/'

        output_filenames = mapping_to_ref.set_ref_output_filenames(self.sample, self.ref_name, self.tmp_folder, self.out_folder)

        self.left_sorted = output_filenames['left_sorted']
        self.right_sorted = output_filenames['right_sorted']
        self.filenames = output_filenames

    def test_create_bed_files_01(self):

        test_intersect_file, test_closest_file, test_left_up, test_right_up = mapping_to_ref.create_bed_files(self.left_sorted, self.right_sorted, self.filenames, 6, '100')

        self.assertTrue(os.path.isfile(test_intersect_file))
        self.assertTrue(os.path.isfile(test_closest_file))
        self.assertTrue(os.path.isfile(test_left_up))
        self.assertTrue(os.path.isfile(test_right_up))
