import luigi

class runcrtm(luigi.Task):
    def output(self):
        return luigi.LocalTarget('profiles/profile2d4_2019_dorain_gfs_output_test_snow_0.nc')
    def run(self):
        print("hello to run crtm")

class Plot_aero(luigi.Task):

    def output(self):
        return luigi.LocalTarget('aero.png')

    def run(self):
        import plot_pdf
        f="profiles/profile2d4_2019_dorain_gfs_output_test_aero_*round3*.nc"
        fig_name="aero.png"
        plot_pdf.plot(f,fig_name)


class Plot_angle0(luigi.Task):
    def requires(self):
        return runcrtm()

    def output(self):
        return luigi.LocalTarget('angle_0.png')

    def run(self):
        import plot_pdf
        f="profiles/profile2d4_2019_dorain_gfs_output_test_angle0*.nc"
        fig_name="angle_0.png"
        plot_pdf.plot(f,fig_name)

class Plot_angle1(luigi.Task):
    def requires(self):
        return runcrtm()

    def output(self):
        return luigi.LocalTarget('angle_1.png')

    def run(self):
        import plot_pdf
        f="profiles/profile2d4_2019_dorain_gfs_output_test_angle1*.nc"
        fig_name="angle_1.png"
        plot_pdf.plot(f,fig_name)
        
class Plot_angle2(luigi.Task):
    def requires(self):
        return runcrtm()

    def output(self):
        return luigi.LocalTarget('angle_2.png')

    def run(self):
        import plot_pdf
        f="profiles/profile2d4_2019_dorain_gfs_output_test_angle2*.nc"
        fig_name="angle_2.png"
        plot_pdf.plot(f,fig_name)
