# Imports to support python 3 compatibility

import os.path
import sys

sys.path.insert(
    0,
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "LINZ", "Geodetic"),
)

import Ellipsoid
import fileunittest


class EllipsoidTestCase(fileunittest.TestCase):
    def test001_Params(self):
        """
        Test handling of ellipsoid parameters
        """
        ell = Ellipsoid.Ellipsoid(6378101, 297.23)
        self.check("001: Semi major axis", ell.a)
        self.check("001: Semi minor axis", ell.b)
        self.check("001: Flattening", ell.rf)

    def test002_Convert(self):
        """
        Conversion geodetic <-> XYZ
        """
        ell = Ellipsoid.Ellipsoid(6378101, 297.23)
        xyz = ell.xyz(172.0, -41.0)
        self.check("002: XYZ from lat/lon", xyz)
        xyz = ell.xyz(172.0, -41.0, 238.0)
        self.check("002: XYZ from lat/lon/hgt", xyz)
        xyz = ell.xyz([172.0, -41.0, 238.0])
        self.check("002: XYZ from [lat,lon,hgt] list", xyz)
        llh = ell.geodetic(xyz)
        self.check("002: LLH from XYZ", llh)
        llh = ell.geodetic([xyz, xyz])
        self.check("002: LLH from multiple XYZ", llh)

    def test003_Calcs(self):
        ell = Ellipsoid.Ellipsoid(6378101, 297.23)
        enu = ell.enu_axes(165.0, -23.0)
        self.check("003: enu_axes", enu.tolist())
        mpd = ell.metres_per_degree(-165.0, -23)
        self.check("003: metres_per_degree lat/lon", mpd)
        mpd = ell.metres_per_degree(165.0, -23.0, 1234.0)
        self.check("003: metres_per_degree lat/lon/hgt", mpd)

    def test004_Grs80(self):
        ell = Ellipsoid.GRS80
        self.check("004: GRS80 Semi major axis", ell.a)
        self.check("004: GRS80 Semi minor axis", ell.b)
        self.check("004: GRS80 Flattening", ell.rf)


if __name__ == "__main__":
    fileunittest.main()
