#!/usr/bin/env python
import re
import pathlib

# Grouping reads and classes for storing this information
READ_PAIRING_REGEXS = [
        # Matches:
        #     prefix_R1.fastq prefix_R2.fastq, or
        #     prefix_R1.fastq.gz prefix_R2.fastq.gz
        re.compile(r'^(.+?)_R[12]\.(fastq(?:\.gz)?)$'),

        # Matches:
        #     prefix_1.fastq prefix_2.fastq, or
        #     prefix_1.fastq.gz prefix_2.fastq.gz
        re.compile(r'^(.+?)_[12]\.(fastq(?:\.gz)?)$'),

        # Matches:
        #     prefix.fastq prefix.fastq, or
        #     prefix.fastq.gz prefix.fastq.gz
        # Note: Final attempt, will match single readsets
        re.compile(r'^(.+?)\.(fastq(?:\.gz)?)$')
                      ]

class ReadSet():

    def __init__(self, prefix, suffix, filepath):
        # Set variables
        self.prefix = prefix
        self.suffix = suffix
        self.filepath = filepath

        # Precalculate some useful things
        self.no_ext = self.filepath.name.replace('.%s' % self.suffix, '')

    def __str__(self):
        return str(self.filepath)


class ReadGroup():

    def __init__(self, prefix, unpaired=None, forward=None, reverse=None):

        # Check that we don't have conflicting arguments
        if bool(forward) ^ bool(reverse):
            raise ValueError('You must pass both forward and reverse')

        if unpaired and forward:
            raise ValueError('You cannot pass both unpaired and forward/ reverse reads')

        # Ensure that all reads have the same prefix and, where applicable, the same suffix
        if forward and reverse:
            assert forward.prefix == reverse.prefix
            assert forward.suffix == reverse.suffix

        # Set the instance variables
        self.prefix = prefix
        self.unpaired = unpaired
        self.forward = forward
        self.reverse = reverse


    # Convenience functions
    @property
    def unpaired_fp(self):
        return self.unpaired_reads.filepath


    @property
    def reverse_fp(self):
        return self.reverse_reads.filepath


    @property
    def forward_fp(self):
        return self.forward_reads.filepath


class ReadGroups():

    def __init__(self, paired, unpaired):
        self.paired = paired
        self.unpaired = unpaired


    def all_groups(self):
        yield from self.paired
        yield from self.unpaired


def group_reads(short_read_fps):
    # Get a map of short read prefixes to ReadSet instances
    read_map = create_prefix_map(short_read_fps)

    # With the remaining short reads, create read groups
    paired_reads = list()
    unpaired_reads = list()

    for prefix, read_sets in read_map.items():
        if len(read_sets) == 1:
            read_group = ReadGroup(prefix=prefix, unpaired=read_sets[0])
            unpaired_reads.append(read_group)

        elif len(read_sets) == 2:
            forward_reads, reverse_reads = sorted(read_sets, key=lambda k: k.filepath)
            read_group = ReadGroup(prefix=prefix, forward=forward_reads, reverse=reverse_reads)
            paired_reads.append(read_group)

        else:
            # Something has gone wrong
            msg_str = ('Too many reads with the same prefix, expected '
                        'either two or one read sets but got: %s')
            raise ValueError(msg_str % ', '.join([str(rd.filepath) for rd in read_sets]))


    # Return instance of the ReadGroups namedtuple
    return ReadGroups(paired_reads, unpaired_reads)


def create_prefix_map(read_fps):
    # Using this method should be O(n * m) at worst;
    # n: number of reads; m: number of regexs
    read_map = dict()

    for read_fp in read_fps:
        # Find suitable regex
        for regex in READ_PAIRING_REGEXS:
            # Apply regex
            re_result = regex.match(read_fp.name)

            # If it works, break
            if re_result:
                break
        else:
            # If no valid regex found, exit and report
            raise ValueError('No regex found for %s (%s)' % (read_fp.name, read_fp))

        # Create and add read_attr to read_pairs using the common prefix as a key
        read_set = ReadSet(re_result.group(1), re_result.group(2), read_fp)
        try:
            read_map[re_result.group(1)].append(read_set)
        except KeyError:
            read_map[re_result.group(1)] = [read_set]

    return read_map

# for testing
def main():
    # groups reads into paired and unpaired
    test = group_reads([pathlib.Path("~/Desktop/ismap_v2/reads/9262_1#29_1.fastq.gz"), pathlib.Path("~/Desktop/ismap_v2/reads/9262_1#29_2.fastq.gz")])
    # ISMapper will only want to work with paired reads
    print(test.paired[0].forward.filepath)
    for group in test.paired:
        print(group.prefix)
        print(group.forward)
        print(group.reverse)
    # But should probably output something sensible if the reads are unpaired
    #print(test.unpaired)
    if not test.unpaired:
        print(test.unpaired)
    # this will show us ALL groups
    #print(test.all_groups())

    # There were no unpaired examples in the first set, here is an example
    #test2 = group_reads([pathlib.Path("~/Desktop/ismap_v2/reads/9262_1#29_1.fastq.gz")])
    #print(test2.unpaired)
    #for group in test2.unpaired:
    #    print(group.prefix)
    #    print(group.unpaired)


if __name__ == '__main__':
    main()