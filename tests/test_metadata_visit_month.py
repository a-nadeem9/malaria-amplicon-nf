import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'runtime' / 'src'))
from simplseq.metadata import enrich_rows_with_metadata


class NominalVisitMonthTests(unittest.TestCase):
    def test_scheduled_month_and_actual_date_are_both_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'metadata.csv'
            path.write_text('participant_id,month,collection_date,pcr\nA,March,2022-04-01,Pos\nA,April,2022-04-29,Neg\n')
            rows=[{'sample_id':'A_March_Rep1','participant_id':'A','collection_date':'2022-03-27','collection_date_source':'inferred_month'},
                  {'sample_id':'A_April_Rep1','participant_id':'A','collection_date':'2022-04-27','collection_date_source':'inferred_month'}]
            out=enrich_rows_with_metadata(rows,path)
            self.assertEqual([r['collection_date'] for r in out],['2022-04-01','2022-04-29'])
            self.assertEqual([r['metadata_pcr'] for r in out],['Pos','Neg'])

    def test_exact_date_still_takes_precedence_and_conflicts_are_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'metadata.csv'
            path.write_text('participant_id,month,collection_date,pcr\nA,March,2022-04-01,Pos\nA,April,2022-04-29,Neg\n')
            rows=[{'sample_id':'A_March','participant_id':'A','collection_date':d,'collection_date_source':'filename_exact_date'}
                  for d in ['2022-04-29','2022-04-02']]
            out=enrich_rows_with_metadata(rows,path)
            self.assertEqual(out[0]['metadata_pcr'],'Neg')
            self.assertEqual(out[1]['metadata_match_status'],'date_conflict')
            self.assertEqual(out[1]['collection_date'],'2022-04-02')

    def test_conflicting_annotations_keep_only_unanimous_date_and_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'metadata.csv'
            path.write_text('participant_id,month,collection_date,pcr,status\nA,May,2022-05-25,Pos,Planned\nA,May,2022-05-25,Neg,Clinical\n')
            out=enrich_rows_with_metadata([{'sample_id':'A_May','participant_id':'A','collection_date':'2022-05-27','collection_date_source':'inferred_month'}],path)[0]
            self.assertEqual(out['collection_date'],'2022-05-25')
            self.assertEqual(out['metadata_match_status'],'ambiguous')
            self.assertEqual(out['metadata_pcr'],'')
            self.assertEqual(out['metadata_status'],'')

if __name__=='__main__':unittest.main()
