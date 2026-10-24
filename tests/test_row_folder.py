from tcspv3.__main__ import row_folder_name

def test_row_folder_name():
    assert row_folder_name(2, 'SrTiO3') == 'row_2_SrTiO3'
    assert row_folder_name(3, 'Fe2 O3') == 'row_3_Fe2_O3'
    assert row_folder_name(4, '../') == 'row_4_empty'
    assert row_folder_name(5, '') == 'row_5_empty'


def test_clear_previous_removes_only_prediction_files(tmp_path):
    from tcspv3.__main__ import clear_previous
    for name in ['prediction.json', 'oxidation.json', 'tcspv3_001.cif', 'tcspv3_020.cif', 'native/tcspv3_001.cif', 'notes.txt']:
        (tmp_path/name).parent.mkdir(exist_ok=True)
        (tmp_path/name).write_text('x')
    clear_previous(tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == ['notes.txt']
