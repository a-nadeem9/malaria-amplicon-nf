"""Boundary cases for the manuscript's observed first/last formula."""
import sys
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'runtime' / 'src'))
from simplseq.clonal import calculate_clonal_dynamics


class ClonalDynamicsTests(unittest.TestCase):
    def calculate(self, rows):
        return calculate_clonal_dynamics(pd.DataFrame(rows, columns=['participant_id','collection_date','locus','allele']))

    def test_first_last_not_reappearance_and_replicates(self):
        rows=[('A','2022-02-01','AMA1','a'),('A','2022-02-01','AMA1','a'),
              ('A','2022-03-01','AMA1','b'),('A','2022-04-01','AMA1','a')]
        tables, summary=self.calculate(rows)
        row=tables['turnover_metrics'].iloc[0]
        self.assertEqual((row.acquisitions,row.clearances,row.genotyped_intervals),(1,1,2))
        self.assertEqual(row.turnover_per_interval,1)
        self.assertEqual(tables['adjacent_visit_transitions'].acquired.sum(),2)
        self.assertEqual(summary['eligible_participant_locus_rows'],1)

    def test_missing_and_single_visit_are_not_absences(self):
        tables, summary=self.calculate([('A','2022-02-01','AMA1','a'),('A','','AMA1','b')])
        self.assertEqual(summary['excluded_calls'],1)
        self.assertEqual(summary['eligible_participant_locus_rows'],0)
        self.assertIsNone(summary['median_turnover_per_interval'])
        self.assertEqual(tables['turnover_metrics'].iloc[0].turnover_per_interval,0)

    def test_other_locus_does_not_change_observation_bounds(self):
        tables, _=self.calculate([('A','2022-02-01','AMA1','a'),('A','2022-04-01','AMA1','a'),
                                  ('A','2022-03-01','CSP','a')])
        rows=tables['turnover_metrics'].set_index('locus')
        self.assertEqual(rows.loc['AMA1','genotyped_intervals'],1)
        self.assertEqual(rows.loc['AMA1','turnover_per_interval'],0)
        self.assertFalse(rows.loc['CSP','eligible_longitudinal'])

if __name__=='__main__':unittest.main()
