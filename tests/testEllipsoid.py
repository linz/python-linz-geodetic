# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
import os.path
sys.path.insert(0,os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'LINZ','Geodetic'))

print(sys.path[0])
from LINZ import fileunittest
import Ellipsoid


class EllipsoidTestCase( fileunittest.TestCase ):

    def test001_Params( self ):
        '''
        Test handling of ellipsoid parameters
        '''
        ell=Ellipsoid.Ellipsoid(6378101,297.23)
        self.check('Semi major axis',ell.a)
        self.check('Semi minor axis',ell.b)
        self.check('Flattening',ell.rf)

    def test002_Params( self ):
        '''
        Conversion geodetic <-> XYZ
        '''
        ell=Ellipsoid.Ellipsoid(6378101,297.23)
        xyz=ell.xyz(172.0,-41.0)
        self.check('XYZ from lat/lon',xyz)
        xyz=ell.xyz(172.0,-41.0,238.0)
        self.check('XYZ from lat/lon/hgt',xyz)
        xyz=ell.xyz([172.0,-41.0,238.0])
        self.check('XYZ from [lat,lon,hgt] list',xyz)

if __name__ == "__main__":
    fileunittest.main()


