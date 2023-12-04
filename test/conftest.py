import pytest

@pytest.fixture
def test_obsids():
    obsids = ['0411980401', # QPE Chakraborty 2021
              '0891800601', # QPE Chakraborty 2021
              '0112570701', # M31 Type 1a X-ray burst Pietsch 2005
              '0001730201', # Provided in examples
              '0002970201', # Provided in examples
              '0831790701'  # Provided in examples
             ]
    return obsids

