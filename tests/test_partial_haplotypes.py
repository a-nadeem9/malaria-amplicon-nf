"""Unknown sequence must not become SNP evidence or full-locus identity."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RT = ROOT / 'runtime'
sys.path.insert(0, str(RT / 'src'))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


cigar = module('partial_cigar', RT / 'workflow/scripts/ASV_to_CIGAR.py')
analysis = module('partial_analysis', RT / 'workflow/scripts/prepare_analysis_table.py')


class PartialComparisonTests(unittest.TestCase):
    def parse(self, root, records, **kwargs):
        path = root / 'toy.msa'
        path.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in records))
        coverage = {}
        result = cigar.parse_alignment(path, amplicon='Toy', min_homopolymer_length=0,
                                       coverage=coverage, **kwargs)
        return result, coverage

    def test_same_region_prevents_wildcard_false_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            records = [('Toy','ACGTACGT'), ('full_ref','ACGTACGT'),
                       ('visible_variant','ACATACGT'), ('partial_variant','ACNNACGT')]
            calls, coverage = self.parse(root, records)
            # The ambiguous site is excluded for EVERY allele, explicitly
            # defining a reduced locus rather than imputing its missing base.
            self.assertEqual(len(set(calls.values())), 1)
            self.assertTrue(calls['full_ref'].startswith('partial@'))
            self.assertEqual(coverage['Toy']['callable_intervals'], [[1,2],[5,8]])
            self.assertEqual(coverage['Toy']['asvs']['visible_variant']['raw_cigar'], '3A')
            reversed_calls, reversed_coverage = self.parse(root, records[:1]+list(reversed(records[1:])))
            self.assertEqual(calls, reversed_calls)
            self.assertEqual(coverage['Toy']['comparison_id'], reversed_coverage['Toy']['comparison_id'])

    def test_known_differences_and_full_observations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            full, coverage = self.parse(root, [('Toy','ACGTACGT'),('ref','ACGTACGT'),('snp','ATGTACGT')])
            self.assertEqual(full, {'ref':'.','snp':'2T'})
            self.assertEqual(coverage, {})
            partial, _ = self.parse(root, [('Toy','ACGTACGT'),('ref','ACNNACGT'),('snp','ATNNACGT')])
            self.assertTrue(partial['ref'].endswith(':.'))
            self.assertTrue(partial['snp'].endswith(':2T'))

    def test_no_common_information_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, 'No shared callable'):
                self.parse(Path(tmp), [('Toy','ACGTACGT'),('left','ACGTNNNN'),('right','NNNNACGT')])

    def test_boundary_deletion_is_not_called_as_a_shorter_event(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls, coverage = self.parse(Path(tmp), [('Toy','ACGTACGT'),('missing','ACGNNNGT'),('deletion','AC---CGT')])
            self.assertTrue(all(call.endswith(':.') for call in calls.values()))
            self.assertTrue(coverage['Toy']['boundary_event_mask_driver_asvs'])

    def test_mixed_regions_and_full_partial_comparisons_are_rejected(self):
        one, two = '1'*32, '2'*32
        for cols in ([f'CSP,partial@{one}:.',f'CSP,partial@{two}:31A'],
                     [f'CSP,partial@{one}:.','CSP,.']):
            with self.assertRaisesRegex(ValueError, 'different observation regions'):
                analysis.parse_allele_columns(cols)
        parsed = analysis.parse_allele_columns([f'CSP,partial@{one}:.',f'CSP,partial@{one}:31A','AMA,.'])
        self.assertEqual(len(parsed),3)


@unittest.skipUnless(shutil.which('Rscript') and shutil.which('muscle'), 'Workflow R and MUSCLE are required')
class PartialPipelineTests(unittest.TestCase):
    def run_command(self, *args):
        result = subprocess.run(list(map(str,args)), capture_output=True, text=True, timeout=120,
                                env={**os.environ,'SIMPLSEQ_THREADS':'1'})
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_standalone_r_bridges_enforce_the_same_comparison_region(self):
        one, two = '1'*32, '2'*32
        cases = [([f'CSP,partial@{one}:.', f'CSP,partial@{one}:31A', 'AMA,.'], None),
                 ([f'CSP,partial@{one}:.', f'CSP,partial@{two}:31A'], 'different observation regions'),
                 ([f'CSP,partial@{one}:.', 'CSP,.'], 'different observation regions'),
                 (['CSP,partial@invalid:.'], 'Malformed partial-region')]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'samples.csv').write_text('sample_id,participant_id,collection_date,replicate\nsample1,person1,2026-09-10,1\n')
            for script in ('simplseq_to_dinemites.R', 'simplseq_to_dcifer.R'):
                for columns, error in cases:
                    with self.subTest(script=script, columns=columns):
                        (root/'counts.tsv').write_text('sample\t'+'\t'.join(columns)+'\nsample1\t'+'\t'.join('1000' for _ in columns)+'\n')
                        result = subprocess.run(['Rscript', str(RT/'workflow/scripts'/script),
                            '--cigar', str(root/'counts.tsv'), '--samples', str(root/'samples.csv'),
                            '--out', str(root/'prepared.tsv')], cwd=root, capture_output=True, text=True, timeout=60)
                        if error:
                            self.assertNotEqual(result.returncode, 0)
                            self.assertIn(error, result.stdout+result.stderr)
                        else:
                            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                            self.assertIn(f'partial@{one}:', (root/'prepared.tsv').read_text())

    def test_adjustment_mapping_and_cigar_with_real_reference_controls(self):
        with (RT/'reference/amplicons_noprimers.fasta').open() as handle:
            refs={r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
        cases={}
        for locus, left, right in (('CSP',125,130),('TRAP',130,130)):
            ref=refs[locus]
            other=next(base for base in 'ACGT' if base!=ref[30])
            edge=next(base for base in 'ACGT' if base!=ref[122])
            variants={'ref':ref,'snp':ref[:30]+other+ref[31:],
                      'ins':ref[:30]+other+ref[30:], 'del':ref[:30]+ref[31:],
                      'edge':ref[:122]+edge+ref[123:]}
            for name, seq in variants.items():
                cases[f'{locus}_{name}']=seq[:left]+'N'*10+seq[-right:]
            cases[f'{locus}_short']=ref[:left-10]+'N'*10+ref[-right+10:]
            cases[f'{locus}_edge_short']=variants['edge'][:left-10]+'N'*10+variants['edge'][-right+10:]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            columns=list(dict.fromkeys(cases.values()))
            with (root/'input.tsv').open('w') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['sample',*columns])
                for case,seq in cases.items():w.writerow([case,*[1000 if s==seq else 0 for s in columns]])
            with (root/'samples.csv').open('w') as f:
                w=csv.writer(f);w.writerow(['sample_id','biological_sample_id','sample_type','participant_id','collection_date'])
                for case in cases:w.writerow([case,case,'sample',case,'2026-09-10'])
            self.run_command('Rscript',RT/'workflow/scripts/adjustASV.R','-s',root/'input.tsv',
                             '-ref',RT/'reference/amplicons_noprimers.fasta','-o',root/'corrected.tsv')
            self.run_command('Rscript',RT/'workflow/scripts/postProc_dada2.R','-s',root/'seqtab_corrected.tsv',
                             '--samples',root/'samples.csv','-ref',RT/'reference/amplicons_noprimers.fasta',
                             '--snv_filter',RT/'reference/snv_filters.txt','-o',root/'mapped.tsv','--fasta')
            with (root/'mapped.tsv').open() as f: mapped={r['hapid']:r for r in csv.DictReader(f,delimiter='\t')}
            with (root/'ASVSeqs.fasta').open() as handle:
                sequences={str(r.seq):r.id for r in SeqIO.parse(handle,'fasta')}
            with (root/'corrected.tsv').open() as f: corrected={r['ASV']:r['correctedASV'] for r in csv.DictReader(f,delimiter='\t')}
            self.run_command(sys.executable,RT/'workflow/scripts/ASV_to_CIGAR.py',root/'ASVSeqs.fasta',root/'mapped.tsv',
                             root/'seqtab_corrected.tsv',root/'counts.tsv','--asv_to_cigar',root/'cigars.tsv',
                             '--amp_db',RT/'reference/amplicons_noprimers.fasta','--alignments',root/'alignments','-p','5')
            with (root/'cigars.tsv').open() as f: calls={r['ASV']:r['CIGAR'] for r in csv.DictReader(f,delimiter='\t')}
            for locus in ('CSP','TRAP'):
                ids={name:sequences[corrected[cases[f'{locus}_{name}']]] for name in ('ref','snp','ins','del','edge','short','edge_short')}
                self.assertEqual(mapped[ids['ref']]['snv_dist_from_3D7'],'0')
                self.assertEqual(mapped[ids['snp']]['snv_dist_from_3D7'],'1')
                self.assertTrue(all(int(row['snv_dist_from_3D7'])>=0 for row in mapped.values()))
                self.assertTrue(all(row['mapping_status']=='MAPPED' for row in mapped.values()))
                self.assertTrue(all(row['strain']=='N' for row in mapped.values()))
                self.assertEqual(calls[ids['ref']],calls[ids['short']])
                self.assertEqual(calls[ids['ref']],calls[ids['edge']])
                self.assertEqual(calls[ids['edge']],calls[ids['edge_short']])
                self.assertNotEqual(calls[ids['ref']],calls[ids['snp']])
                self.assertRegex(calls[ids['ins']],r':31I=[ACGT]+$')
                self.assertRegex(calls[ids['del']],r':31D=[ACGT]+$')
            with (root/'counts.tsv').open() as f:
                rows=list(csv.reader(f,delimiter='\t'))
            self.assertEqual(sum(int(n) for row in rows[1:] for n in row[1:]),14000)
            coverage=json.loads((root/'counts.tsv.coverage.json').read_text())
            self.assertEqual(set(coverage['partial_loci']),{'CSP','TRAP'})


if __name__=='__main__': unittest.main()
