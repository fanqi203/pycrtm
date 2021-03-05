import plot_pdf
import os.path
from os import path

def test_answer():
    f="profiles/profile2d4_2019_dorain_gfs_output_test_angle0*.nc"
    fig_name="angle_0_tmp.png"
    if path.exists(fig_name):
        os.remove(fig_name)
    plot_pdf.plot(f,fig_name)
    assert path.exists(fig_name)


    
