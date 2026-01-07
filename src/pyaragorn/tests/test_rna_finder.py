import itertools
import re
import unittest

from .. import RNAFinder, TMRNAGene, TRNAGene
from . import data


_TRNA_RX = re.compile(r"^(\d+)\s+tRNA-([A-Za-z]{3})\s+(c?)\[(\d+),(\d+)\]\s+([\d.]+)\s+(\d+)\s+\(([a-z]{2,4})\)")
_TMRNA_RX = re.compile(r"^(\d+)\s+tmRNA\s+(c?)\[(\d+),(\d+)]\s+([\d.]+)\s+(\d+),(\d+)\s+([A-Z\*]+)")

def batched(iterable, n, *, strict=False):
    # batched('ABCDEFG', 3) → ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    iterator = iter(iterable)
    for batch in iter(lambda: tuple(itertools.islice(iterator, n)), ()):
        if strict and len(batch) != n:
            raise ValueError('batched(): incomplete batch')
        yield batch

class TestRNAFinder(unittest.TestCase):
    
    def test_default(self):
        record = data.load_record("CP001621.fna.gz")
        lines = data.load_text("CP001621.default.txt").splitlines()
        
        finder = RNAFinder(translation_table=11)
        genes = finder.find_rna(str(record.seq))

        for gene, expected in itertools.zip_longest(genes, batched(lines[2:], 3)):
            self.assertIsNotNone(gene)
            self.assertIsNotNone(expected)
            result, seq, ss = expected
            if gene.type == "tRNA":
                matched = _TRNA_RX.match(result)
                _, aa, complement, begin, end, energy, offset, anticodon = matched.groups()
                self.assertEqual(gene.amino_acid, aa)
                self.assertEqual(gene.begin, int(begin))
                self.assertEqual(gene.end, int(end))
                self.assertEqual(gene.anticodon_offset, int(offset))
                self.assertEqual(gene.anticodon_length, len(anticodon))
                self.assertEqual(gene.anticodon, anticodon)
                self.assertEqual(gene.strand, -1 if complement == "c" else +1)
                self.assertAlmostEqual(gene.energy, float(energy), places=1)
                self.assertEqual(gene.sequence().lower(), seq)
            elif gene.type == "tmRNA":
                matched = _TMRNA_RX.match(result)
                _, complement, begin, end, energy, orf_start, orf_end, peptide = matched.groups()
                self.assertEqual(gene.begin, int(begin))
                self.assertEqual(gene.end, int(end))
                self.assertEqual(gene.orf_offset, int(orf_start))
                self.assertEqual(gene.orf_offset + gene.orf_length, int(orf_end))
                self.assertEqual(gene.peptide(), peptide)
                self.assertEqual(gene.strand, -1 if complement == "c" else +1)
                self.assertAlmostEqual(gene.energy, float(energy), places=1)
                # self.assertEqual(gene.sequence().lower(), seq) # TODO

    def test_trna(self):
        record = data.load_record("CP001621.fna.gz")
        finder = RNAFinder(translation_table=11, tmrna=False, trna=True)
        for gene in finder.find_rna(str(record.seq)):
            self.assertIsInstance(gene, TRNAGene)
    
    def test_tmrna(self):
        record = data.load_record("CP001621.fna.gz")
        finder = RNAFinder(translation_table=11, tmrna=True, trna=False)
        for gene in finder.find_rna(str(record.seq)):
            self.assertIsInstance(gene, TMRNAGene)

    def test_ps95(self):
        record = data.load_record("CP001621.fna.gz")
        lines = data.load_text("CP001621.ps95.txt").splitlines()
        
        finder = RNAFinder(translation_table=11, threshold_scale=0.95)
        genes = finder.find_rna(str(record.seq))

        for gene, expected in itertools.zip_longest(genes, batched(lines[2:], 3)):
            self.assertIsNotNone(gene)
            self.assertIsNotNone(expected)
            result, seq, ss = expected
            if gene.type == "tRNA":
                matched = _TRNA_RX.match(result)
                _, aa, complement, begin, end, energy, offset, anticodon = matched.groups()
                self.assertEqual(gene.amino_acid, aa)
                self.assertEqual(gene.begin, int(begin))
                self.assertEqual(gene.end, int(end))
                self.assertEqual(gene.anticodon_offset, int(offset))
                self.assertEqual(gene.anticodon_length, len(anticodon))
                self.assertEqual(gene.anticodon, anticodon)
                self.assertEqual(gene.strand, -1 if complement == "c" else +1)
                self.assertAlmostEqual(gene.energy, float(energy), places=1)
                self.assertEqual(gene.sequence().lower(), seq)
            elif gene.type == "tmRNA":
                matched = _TMRNA_RX.match(result)
                _, complement, begin, end, energy, orf_start, orf_end, peptide = matched.groups()
                self.assertEqual(gene.begin, int(begin))
                self.assertEqual(gene.end, int(end))
                self.assertEqual(gene.orf_offset, int(orf_start))
                self.assertEqual(gene.orf_offset + gene.orf_length, int(orf_end))
                self.assertEqual(gene.peptide(), peptide)
                self.assertEqual(gene.strand, -1 if complement == "c" else +1)
                self.assertAlmostEqual(gene.energy, float(energy), places=1)
                # self.assertEqual(gene.sequence().lower(), seq) # TODO