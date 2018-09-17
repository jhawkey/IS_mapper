import unittest
import pathlib
import ismap_v2
import read_grouping
import run_commands
import mapping_to_ref
import mapping_to_query
import os
import create_output
from Bio import SeqIO
import filecmp
import shutil

#TODO: Provide empty bed files (try all combinations) and ensure error messages are correct, temp files deleted, program exits nicely
#TODO: provide misformatted gbk ref (LOCUS too long, qualifiers wrong) -> error messages, tmp dir, program exits
#TODO: provide misformatted fasta query or queries
#TODO: provide read names which can't be paired
#TODO: provide read names formatted differently to check pairing

#TODO: run all 8 Ab samples with each ref gbk, try a run where both gbks provided as either multi-genbank or two entries

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

        test_left, test_right, is_output_folder, tmp_output_folder = ismap_v2.map_to_is_query(self.sample, self.is_query, self.output)

        # verify that the read files are actually being made
        self.assertTrue(os.path.isfile(test_left))
        self.assertTrue(os.path.isfile(test_right))

        # verify that the output folders are correct
        self.assertTrue(is_output_folder == os.path.join(self.output, 'ISAba1'))
        self.assertTrue(tmp_output_folder == os.path.join(self.output, 'ISAba1', 'tmp'))

        # verify that we have the same number of reads, and the reads are what we expect using filecmp
        # increase the buffersize so we check more than the first 8Kb of the file, check first 8Mb
        filecmp.BUFSIZE = 1024 * 10
        # set shallow to be false so that it actually checks the contents of the file, not just os.stat
        gold_left = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_ISAba1_left_final.fastq'
        gold_right = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_ISAba1_right_final.fastq'
        self.assertTrue(filecmp.cmp(test_left, gold_left ,shallow=False))
        self.assertTrue(filecmp.cmp(test_right, gold_right, shallow=False))

    def tearDown(self):

        # remove directory containing files after test
        shutil.rmtree(self.output)

class TestEmptyFiles(unittest.TestCase):

    def setUp(self):
        self.left_flanking = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_ISAba1_left_final.fastq'
        self.empty_right_flanking = '/Users/jane/Desktop/ismap_v2/gold_standard_files/empty_files/9262_1#29_ISAba1_right_final.fastq'

        ref_seqs = ismap_v2.get_sequences(['/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk'], 'genbank')
        self.ref = ref_seqs[0]

        self.query_single = ['/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta']
        self.query_records = ismap_v2.get_sequences(self.query_single, 'fasta')
        self.is_query = self.query_records[0]

        self.sample = '9262_1#29'
        self.main_out_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/'
        self.tmp_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/tmp'
        self.out_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/'
        # make the output directories
        if not os.path.exists(self.tmp_folder):
            os.makedirs(self.tmp_folder)

    def test_empty_read_file(self):
        # providing one empty read file from the query stage
        # test should just run and return no errors if everything works fine
        test_filenames = mapping_to_ref.map_to_ref_seq(self.ref, self.sample, self.left_flanking,
                                                       self.empty_right_flanking, self.tmp_folder, self.out_folder, '1', False)

        test_filenames = mapping_to_ref.create_bed_files(test_filenames, 6, '100')

        create_output.create_typing_output(test_filenames, self.ref, self.is_query, 0.9, 1.1, self.tmp_folder, self.sample)

    def tearDown(self):
        # remove directory containing files after test
        shutil.rmtree(self.main_out_folder)

class TestRefMapping(unittest.TestCase):

    def setUp(self):

         self.left_flanking = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_ISAba1_left_final.fastq'
         self.right_flanking = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_ISAba1_right_final.fastq'

         ref_seqs = ismap_v2.get_sequences(['/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk'], 'genbank')
         self.ref = ref_seqs[0]
         self.sample = '9262_1#29'
         self.main_out_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/'
         self.tmp_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/tmp'
         self.out_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/'
         # make the output directories
         if not os.path.exists(self.tmp_folder):
            os.makedirs(self.tmp_folder)

    def test_map_to_ref_seq_01(self):

        # test that after mapping to ref the reads in the sorted bam files are correct
        test_filenames = mapping_to_ref.map_to_ref_seq(self.ref, self.sample, self.left_flanking,
                                                       self.right_flanking, self.tmp_folder, self.out_folder, '1', False)
        test_left_sorted_bam = test_filenames['left_sorted']
        test_right_sorted_bam = test_filenames['right_sorted']

        # verify that the resulting bam files are the same as the gold standard using filecmp
        # increase the buffersize so we check more than the first 8Kb of the file, check first 8Mb
        filecmp.BUFSIZE = 1024 * 10

        gold_left_sorted_sam = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_left_CP010781.1.sorted.sam'
        gold_right_sorted_sam = '/Users/jane/Desktop/ismap_v2/gold_standard_files/9262_1#29_right_CP010781.1.sorted.sam'

        # we need to convert the BAM files to SAM files to do the check, as the BAM files contain a header
        # specifying the version of samtools
        # set up samtools
        samtools_runner = mapping_to_query.RunSamtools()
        # run samtools view to convert to SAM
        run_commands.run_command(samtools_runner.view_bam_to_sam(test_left_sorted_bam, test_left_sorted_bam + '.sam'), shell=True)
        run_commands.run_command(samtools_runner.view_bam_to_sam(test_right_sorted_bam, test_right_sorted_bam + '.sam'), shell=True)

        # check if the SAM files match
        self.assertTrue(filecmp.cmp(test_left_sorted_bam + '.sam', gold_left_sorted_sam, shallow=False))
        self.assertTrue(filecmp.cmp(test_right_sorted_bam + '.sam', gold_right_sorted_sam, shallow=False))

    def tearDown(self):
        # remove directory containing files after test
        shutil.rmtree(self.main_out_folder)

class TestCreateBedFiles(unittest.TestCase):

    def setUp(self):
        pass

    def test_create_bed_files_01(self):
        pass

    def tearDown(self):
        pass


class TestSetOutputFilenames(unittest.TestCase):
    def setUp(self):
        self.tmp_folder = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/tmp'
        self.prefix = '9262_1#29'
        self.query = 'ISAba1'
        self.out_dir = '/Users/jane/Desktop/ismap_v2/test_results/9262_1#29/ISAba1/'

    def test_set_output_filenames(self):
        pass

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

class TestCreate_typing_output(unittest.TestCase):


    def setUp(self):
        #self.filenames = {'intersect':'/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP010781.1_intersect.bed',
        #                  'closest': '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP010781.1_closest.bed'}
        #self.ref_gbk_obj = SeqIO.read('/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk', 'genbank')
        self.filenames = {'intersect':'/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP001921.1_intersect.bed',
                          'closest': '/Users/jane/Desktop/ismap_v2/9262_1#29/ISAba1/9262_1#29_CP001921.1_closest.bed'}
        self.ref_gbk_obj = SeqIO.read('/Users/jane/Desktop/ismap_v2/refs/CP001921.gbk', 'genbank')
        self.is_query_obj = SeqIO.read('/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta', 'fasta')
        self.min_range = 0.9
        self.max_range = 1.1

    def test_create_typing_output(self):
        hit_list = create_output.create_typing_output(self.filenames, self.ref_gbk_obj, self.is_query_obj, self.min_range, self.max_range)