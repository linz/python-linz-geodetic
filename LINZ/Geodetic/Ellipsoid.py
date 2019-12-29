#/usr/bin/python3

import numpy as np
import math

class Ellipsoid( object ):

    convergence=1.0e-10

    @staticmethod
    def _cossin( angle ):
        angle=np.radians(angle)
        return np.cos(angle),np.sin(angle)

    @staticmethod
    def enu_axes( lon, lat ):
        '''
        Returns an array defining the east, north, and up unit vectors
        at a specified latitude and longitude

        To convert an xyz offset to an enu offset, use as an example
        
           enu_axes=GRS80.enu_axes(lon,lat)
           denu=enu_axes.dot(dxyz)
           dxyz=enu_axes.T.dot(denu)

        '''
        cln,sln = Ellipsoid._cossin(lon)
        clt,slt = Ellipsoid._cossin(lat)
        ve=np.array([-sln,cln,0])
        vn=np.array([-cln*slt,-sln*slt,clt])
        vu=np.array([clt*cln,clt*sln,slt])
        return np.vstack((ve,vn,vu))

    def __init__( self, a, rf ):
        '''
        Initiallize an ellipsoid based on semi major axis and inverse flattening
        '''
        self._setParams(a,rf)

    @property
    def a( self ): return self._a

    @a.setter
    def a( self, a ): self._setParams( a, self._rf )

    @property
    def rf( self ): return self._rf

    @rf.setter
    def rf( self, rf ): self._setParams( self._a, rf )

    @property
    def b( self ): return self._b

    def _setParams(self,a,rf):
        self._a=float(a)
        self._rf=float(rf)
        self._b=a-a/rf if rf else a
        self._a2=a*a
        self._b2=self._b*self._b
        self._a2b2=self._a2-self._b2

    def xyz( self, lon, lat=None, hgt=None ):
        '''
        Calculate the geocentric X,Y,Z coordinates at a longitude 
        and latitude

        Input is one of 
           lon, lat        Single values or lists of values
           lon, lat, hgt 
           llh             Array of [lon,lat,hgt]
        '''
        single=True
        if lat is None:
            if not isinstance(lon,np.ndarray):
                lon=np.array(lon)
            single=len(lon.shape)==1
            if single:
                lon=lon.reshape((1,lon.size))
            lat=lon[:,1]
            hgt=lon[:,2] if lon.shape[1] > 2 else 0
            lon=lon[:,0]
        if hgt is None:
            hgt=0

        cln,sln = Ellipsoid._cossin(lon)
        clt,slt = Ellipsoid._cossin(lat)
        bsac=np.hypot(self._b*slt,self._a*clt)
        p = self._a2*clt/bsac + hgt*clt
        
        xyz=[p*cln,p*sln,self._b2*slt/bsac+hgt*slt]
        xyz=np.vstack(xyz).transpose()
        if single:
            xyz=xyz[0]
        return xyz

    def metres_per_degree( self, lon, lat, hgt=0 ):
        '''
        Calculate the number of metres per degree east and
        north
        '''
        clt,slt = Ellipsoid._cossin(lat)
        bsac=np.hypot(self._b*slt,self._a*clt)
        p = self._a2*clt/bsac + hgt*clt
        dedln=np.radians(p)
        dndlt=np.radians((self._a2*self._b2)/(bsac*bsac*bsac)+hgt)
        return dedln,dndlt

    def geodetic( self, xyz ):
        '''
        Calculate the longitude, latitude, and height corresponding 
        to a geocentric XYZ coordinate

        Input is one of 
           xyz             Single [x,y,z]
           xyz             Array of [x,y,z]
        '''
        if not isinstance(xyz,np.ndarray):
            xyz=np.array(xyz)
        single=len(xyz.shape)==1
        if single:
            xyz=xyz.reshape((1,xyz.size))
        x,y,z = xyz[:,0],xyz[:,1],xyz[:,2]
        ln=np.arctan2(y,x)
        p=np.hypot(x,y)
        lt=np.arctan2(self._a2*z,self._b2*p)
        for i in range(10):
            lt0=lt
            slt=np.sin(lt)
            clt=np.cos(lt)
            bsac=np.hypot(self._b*slt,self._a*clt)
            lt=np.arctan2(z+slt*self._a2b2/bsac,p)
            if np.all(abs(lt-lt0) < self.convergence):
                break
        h=p*clt+z*slt-bsac
        result=np.degrees(ln),np.degrees(lt),h
        result=np.vstack(result).transpose()
        return result[0] if single else result

    def radii_of_curvature( self, latitude ):
        '''
        Returns the radii of curvature along the meridional and prime
        vertical normal sections.
        '''
        clt,slt = Ellipsoid._cossin(latitude)
        den=math.sqrt(self._a2*clt*clt+self._b2*slt*slt)
        rm=self._a2*self._b2/(den*den*den)
        rn=self._a2/den
        return rm, rn


GRS80 = Ellipsoid(6378137.0,298.257222101)

def main():
    import sys
    import argparse
    import re
    parser=argparse.ArgumentParser(description='Convert Cartesian coordinates <=> Geodetic coordinates')
    parser.add_argument('-e','--ellipsoid',type=float,nargs=2,metavar='A RF',help='Ellipsoid semi-major axis and inverse flattening (default GRS80)')
    parser.add_argument('-x','--xyz',nargs=3,type=float,metavar=('X','Y','Z'),help="XYZ coordinates to transform")
    parser.add_argument('-g','--geodetic',nargs=3,type=float,metavar=('LON','LAT','EHGT'),help="Geodetic coordinates to transform")
    parser.add_argument('-r','--calc_geodetic',action='store_true',help="Calculate geodetic from XYZ (default is calculate XYZ)")
    parser.add_argument('-c','--csv',action='store_true',help="File format CSV - default whitespace delimited")
    parser.add_argument('-zg','--geodetic-column-names',metavar=('LON_COL','LAT_COL','HGT_COL'),nargs=3,help="Column names of X,Y,Z fields - default first three columns")
    parser.add_argument('-zx','--xyz-column-names',metavar=('X_COL','Y_COL','Z_COL'),nargs=3,help="Column names of X,Y,Z fields - default first three columns")
    parser.add_argument('input_file',nargs='?',help="Input file of XYZ coordinates")
    parser.add_argument('output_file',nargs='?',help="Output file of XYZ coordinates")

    args=parser.parse_args()
    input_file=args.input_file
    output_file=args.output_file
    
    if args.xyz is not None and input_file is not None:
        print("Cannot have xyz and input file arguments")
        sys.exit()
    if args.geodetic is not None and input_file is not None:
        print("Cannot have geodetic and input file arguments")
        sys.exit()
    if args.xyz is not None and args.geodetic is not None:
        print("Cannot have xyz and geodetic arguments")

        sys.exit()
    if args.geodetic is None and args.xyz is None and not input_file:
        print("No coordinate input specified - need xyz, geodetic, or input file")
        sys.exit()
    if input_file is not None and output_file is None:
        print("Need the name of an output file")
        sys.exit()

    ell=GRS80
    if args.ellipsoid:
        ell=Ellipsoid(args[0],args[1])

    if args.xyz:
        llh=ell.geodetic(args.xyz)
        print("{0:.9f} {1:.9f} {2:.4f}".format(*llh))
        sys.exit()

    if args.geodetic:
        xyz=ell.xyz(args.geodetic)
        print("{0:.4f} {1:.4f} {2:.4f}".format(*xyz))
        sys.exit()

    tfm = ell.geodetic if args.calc_geodetic else ell.xyz
    gcol=args.geodetic_column_names
    xcol=args.xyz_column_names
    if gcol and not xcol:
        xcol=['X','Y','Z']
    if xcol and not gcol:
        gcol=['Lon','Lat','Hgt']

    incol=xcol if args.calc_geodetic else gcol
    outcol=gcol if args.calc_geodetic else xcol
    colfmt=['{0:.9f}','{0:.9f}','{0:.4f}'] if args.calc_geodetic else ['{0:.4f}']*3

    if args.input_file:
        cols=[0,1,2]
        reqlen=3
        with sys.stdin if input_file=='-' else open(input_file,'r') as fin:
            import csv
            if args.csv:
                freader=csv.reader(fin)
            else:
                def wsreader(f):
                    for l in f:
                        yield l.split()
                freader=wsreader(fin)
            with sys.stdout if output_file=='-' else open(output_file,'w') as fout:
                if args.csv:
                    writerow=csv.writer(fout).writerow
                else:
                    def writerow(row):
                        fout.write('\t'.join(row))
                        fout.write('\n')
                if incol:
                    header=next(freader)
                    header=[x.upper() for x in header]
                    cols=[]
                    for c in incol:
                        c=c.upper()
                        if c not in header:
                            print('Column',c,'is missing from input file header')
                            sys.exit()
                        nc=header.index(c)
                        cols.append(nc)
                        cols.append(header.index(c))
                    reqlen=max(cols)
                    row=list(header)
                    for i,c in zip(cols,outcol):
                        row[i]=c
                    writerow(row)
                for row in freader:
                    if len(row) < reqlen:
                        continue
                    xyz=[float(row[i]) for i in cols]
                    xyzt=tfm(xyz)
                    for c,x,f in zip(cols,xyzt,colfmt):
                        row[c]=f.format(x)
                    writerow(row)

if __name__=="__main__":
    main()
