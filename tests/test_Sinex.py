import sys
import os.path
from datetime import datetime

testdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(testdir), "LINZ"))

import fileunittest
from Geodetic import Sinex


class SinexTestCase(fileunittest.TestCase):

    sinexfile1 = os.path.join(testdir, "test.snx")

    def test001_ReaderSolutionId(self):
        """
        Test basic sinex reading
        """

        snx = Sinex.Reader(SinexTestCase.sinexfile1, covariance=Sinex.COVAR_FULL)
        self.check("test001: solutions", snx.solutions())
        self.check("test001: getSolutionId(KAIK)", snx.getSolutionIds(ptid="KAIK"))
        self.check(
            "test001: getSolutionsIds(KAIK,A)",
            snx.getSolutionIds(ptid="KAIK", ptcode="A"),
        )
        self.check(
            "test001: getSolutionsIds(KAIK,B)",
            snx.getSolutionIds(ptid="KAIK", ptcode="B"),
        )
        # self.check("test001: getSolutionIds(A)",snx.getSolutionIds(ptcode='A'))
        self.check(
            "test001: getSolutionIds(KAIK,1)",
            snx.getSolutionIds(ptid="KAIK", solnid="1"),
        )
        self.check(
            "test001: getSolutionIds(KAIK,2)",
            snx.getSolutionIds(ptid="KAIK", solnid="2"),
        )
        self.check(
            "test001: getSolutionIds(KAIK:A) single", snx.getSolutionIds("KAIK:A")
        )
        self.check(
            "test001: getSolutionIds(KAIK:B) single", snx.getSolutionIds("KAIK:B")
        )
        self.check(
            "test001: getSolutionIds(KAIK:A:1) single", snx.getSolutionIds("KAIK:A:1")
        )
        self.check(
            "test001: getSolutionIds(KAIK:A:2) single", snx.getSolutionIds("KAIK:A:2")
        )

    def test002_ReaderGet(self):
        """
        Test get (no velocities)
        """

        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_FULL, velocities=False
        )
        self.check("test002: get(KAIK)", snx.get("KAIK"))
        self.check("test002: get(solutions(0))", snx.get(snx.solutions()[0]))
        self.check("test002: get(solution(KAIK A))", snx.get("KAIK", ptcode="A"))

    def test003_ReaderGetSolutions(self):
        """
        Test getSolutions (no velocities)
        """
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_FULL, velocities=False
        )
        solutions = snx.solutions()
        self.check("test003: getsolutions(0)", snx.getSolutions(solutions[0:1]))
        self.check("test003: getsolutions(1:3))", snx.getSolutions(solutions[1:3]))
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_STATION, velocities=False
        )
        self.check(
            "test003: getsolutions(0) stn covar", snx.getSolutions(solutions[0:1])
        )
        self.check(
            "test003: getsolutions(1:3)) stn covar", snx.getSolutions(solutions[1:3])
        )
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_NONE, velocities=False
        )
        self.check(
            "test003: getsolutions(0) stn none", snx.getSolutions(solutions[0:1])
        )
        self.check(
            "test003: getsolutions(1:3)) stn none", snx.getSolutions(solutions[1:3])
        )

    def test004_ReaderXYZ(self):
        """
        Test xyz function
        """
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_FULL, velocities=False
        )
        self.check("test004: xyz KAIK", snx.xyz("KAIK"))
        self.check("test004: xyz KAIK tuple", snx.xyz(("KAIK", "A", "1")))
        self.check("test004: xyz KAIK A 1", snx.xyz("KAIK:A:1"))
        self.check("test004: xyz KAIK with covar", snx.xyz("KAIK", covariance=True))
        self.checkRun("test004: xyz ABCD", lambda: snx.xyz("ABCD"))
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_STATION, velocities=False
        )
        self.check("test004: stn covar xyz KAIK", snx.xyz("KAIK"))
        self.check(
            "test004: stn covar xyz KAIK with covar", snx.xyz("KAIK", covariance=True)
        )
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_NONE, velocities=False
        )
        self.check("test004: no covar xyz KAIK", snx.xyz("KAIK"))
        self.checkRun(
            "test004: no covar xyz KAIK with covar",
            lambda: snx.xyz("KAIK", covariance=True),
        )

    def test005_ReaderDXYZ(self):
        """
        Test dxyz function
        """
        snx = Sinex.Reader(
            SinexTestCase.sinexfile1, covariance=Sinex.COVAR_FULL, velocities=False
        )
        self.check("test005: dxyz 1163 KAIK", snx.dxyz("1163", "KAIK"))
        self.check("test005: dxyz KAIK 1163", snx.dxyz("KAIK", "1163"))
        self.checkRun("test005: dxyz KAIK KAIK", lambda: snx.dxyz("KAIK", "KAIK"))
        self.checkRun("test005: dxyz KAIK ABCD", lambda: snx.dxyz("KAIK", "ABCD"))

    def test006_WriteSindex(self):
        filename = os.path.join(testdir, "test006_output.snx")
        with Sinex.Writer(filename) as snxf:
            snxf.addFileInfo(description="Test SINEX.Writer", agency="Land Info NZ")
            snxf.setSolutionStatistics(2.4, 23.8, 180, 6)
            snxf.setObsDateRange(datetime(2010, 2, 5), datetime(2010, 2, 8))
            snxf.addMark("12cd", description="Useful mark for some")
            snxf.addMark("ax1d", monument="Second mark", description=None)
            snxf.addSolution(
                "12cd",
                [-4685480.3689522203, 531054.576640439, -4280819.1694681998],
                [0, 1, 2],
            )
            snxf.addSolution(
                "ax1d",
                [-46854280.3689522203, 531354.576640439, -4280849.1694681998],
                [3, 4, 5],
            )
            snxf.setCovariance(
                [
                    [
                        0.0000003003,
                        -0.0000000264,
                        0.0000002349,
                        0.0000000094,
                        0.0000000012,
                        -0.0000000004,
                        -0.0000000006,
                        0.0000000032,
                    ],
                    [
                        -0.0000000264,
                        0.0000000159,
                        -0.0000000218,
                        -0.0000000006,
                        0.0000000001,
                        0.0000000002,
                        0.000000001,
                        -0.0000000001,
                    ],
                    [
                        0.0000002349,
                        -0.0000000218,
                        0.0000002229,
                        0.0000000062,
                        0.0000000011,
                        -0.0000000016,
                        -0.0000000001,
                        0.0000000027,
                    ],
                    [
                        0.0000000094,
                        -0.0000000006,
                        0.0000000062,
                        0.0000001599,
                        -0.000000014,
                        0.0000001283,
                        -0.000000076,
                        0.0000000078,
                    ],
                    [
                        0.0000000012,
                        0.0000000001,
                        0.0000000011,
                        -0.000000014,
                        0.0000000084,
                        -0.0000000119,
                        0.0000000078,
                        -0.0000000039,
                    ],
                    [
                        -0.0000000004,
                        0.0000000002,
                        -0.0000000016,
                        0.0000001283,
                        -0.0000000119,
                        0.0000001238,
                        -0.0000000612,
                        0.0000000065,
                    ],
                ]
            )
        with open(filename) as sf:
            sinexdata = sf.read()
        sinexdata = sinexdata[:15] + "00:001:12345" + sinexdata[27:]
        self.check("test006: write sinex", sinexdata)
        os.remove(filename)


if __name__ == "__main__":
    fileunittest.main()
