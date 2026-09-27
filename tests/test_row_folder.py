from tcspv3.__main__ import row_folder_name

def test_row_folder_name():
    assert row_folder_name(2, 'SrTiO3') == 'row_2_SrTiO3'
    assert row_folder_name(3, 'Fe2 O3') == 'row_3_Fe2_O3'
    assert row_folder_name(4, '../') == 'row_4_empty'
    assert row_folder_name(5, '') == 'row_5_empty'
