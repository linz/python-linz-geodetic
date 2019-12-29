import sys
import os.path
from datetime import datetime

testdir=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(testdir),'LINZ'))

import fileunittest
from Geodetic import ITRF


class ITRFTestCase( fileunittest.TestCase ):

    def test001_ITRF_parameters( self ):
        '''
        Test ITRF parameters
        '''
        params={p[0]:p[1] for p in ITRF.ITRF_params}
        for r in ('ITRF96','ITRF97','ITRF2000','ITRF2005','ITRF2008','ITRF2014'):
            self.check("test001: {0}".format(r),params.get(r))

    def test002_BWtransform( self ):
        tf=ITRF.Transformation('ITRF96','ITRF2008')
        self.check("test002 BW transform 2000", str(tf))
        tf=ITRF.Transformation('ITRF96','ITRF2008')
        self.check("test002 BW transform 2020",str(tf))
        tf2=tf.atDate(2020.0)
        self.check("test002 BW transform 2020 at 2020",str(tf))
        self.check("test002 BW apply1",tf2.transform([-4761241.541,754106.577,-4162423.201]))
        self.check("test002 BW apply2",tf.transform([-4761241.541,754106.577,-4162423.201],2020.0))
        self.check("test002 BW apply3",tf.transform([-4761241.541,754106.577,-4162423.201]))
        self.check("test002 BW apply4",tf.transformLonLat([171.28,-41.54,34.0]))

if __name__ == "__main__":
    fileunittest.main()


