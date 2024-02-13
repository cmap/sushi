import unittest
from functions.remove_functions import remove_data_based_on_criteria
import pandas as pd
from pandas.testing import assert_frame_equal

class TestRemoveDataBasedOnCriteria(unittest.TestCase):
    def setUp(self):
        # Sample data DataFrame
        self.data_df = pd.DataFrame({
            'project_code': ['EPS001', 'EPS001', 'EPS001', 'EPS002', 'EPS002', 'EPS002'],
            'CCLE_name': ['NCIH358_LUNG', 'SNU840_OVARY', 'SKMEL30_SKIN', 'TE1_OESOPHAGUS', 'TE1_OESOPHAGUS', 'NCIH358_LUNG'],
            'pcr_well': ['A01', 'A02', 'A03', 'A04', 'A04','A05']
        })
        # Sample criteria DataFrame
        self.criteria_df_single = pd.DataFrame({
            'project_code': ['EPS001'],
            'CCLE_name': [pd.NA]
        })
        
        self.criteria_df_multiple = pd.DataFrame({
            'project_code': ['EPS001', 'EPS001'],
            'CCLE_name': ['NCIH358_LUNG', 'SKMEL30_SKIN']
        })
        
        self.criteria_df_bad_column = pd.DataFrame({
            'project_code': ['EPS001','EPS002'],
            'bad_column_header': [1, 2]
        })
        
        self.criteria_df_no_match = pd.DataFrame({
            'project_code': ['MTS024', 'CPS010'],
            'pcr_well': ['B01','B02']
        })
        
        
    def test_remove_single_criteria(self):
        """
        Test data removal based on a single criteria.
        """
        modified_df, n_rows_removed = remove_data_based_on_criteria(self.data_df, self.criteria_df_single)
        expected_df = self.data_df[~self.data_df['project_code'].isin(['EPS001'])]
        
        assert_frame_equal(modified_df.reset_index(drop=True), expected_df.reset_index(drop=True))
        self.assertEqual(n_rows_removed, 3)


    def test_remove_multiple_criteria(self):
        """
        Test data removal based on multiple criteria.
        """
        modified_df, n_rows_removed = remove_data_based_on_criteria(self.data_df, self.criteria_df_multiple)
        expected_df = self.data_df[~((self.data_df['project_code'] == 'EPS001') & (self.data_df['CCLE_name'] == 'NCIH358_LUNG')) &
                                    ~((self.data_df['project_code'] == 'EPS001') & (self.data_df['CCLE_name'] == 'SKMEL30_SKIN'))]
        
        assert_frame_equal(modified_df.reset_index(drop=True), expected_df.reset_index(drop=True))
        self.assertEqual(n_rows_removed, 2)
        
        
    def test_remove_bad_column_headers(self):
        """
        Test data removal with criteria containing bad column headers.
        """
        with self.assertRaises(ValueError) as context:
            _, _ = remove_data_based_on_criteria(self.data_df, self.criteria_df_bad_column)
        
        self.assertTrue('Columns not found in data_df' in str(context.exception))
        
        
    def test_remove_no_matching_data(self):
        """
        Test data removal when criteria match no rows in the data DataFrame.
        """
        with self.assertRaises(ValueError) as context:
            remove_data_based_on_criteria(self.data_df, self.criteria_df_no_match)
        
        self.assertEqual(str(context.exception), "No rows were removed based on the provided criteria.")



if __name__ == '__main__':
    unittest.main(verbosity=2)