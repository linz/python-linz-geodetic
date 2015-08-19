#!/usr/bin/python

import math
import numpy as np
from datetime import date, datetime

# IERS parameters:
#
# Data from http://itrf.ensg.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt
#
# Note : These parameters are derived from those already published in the IERS
# Technical Notes and Annual Reports. The transformation parameters should be
# used with the standard model (1) given below and are valid at the indicated
# epoch.
# 
# 
# : XS :    : X :   : Tx :   :  D   -Rz   Ry : : X :
# :    :    :   :   :    :   :               : :   :
# : YS :  = : Y : + : Ty : + :  Rz   D   -Rx : : Y :                       (1)
# :    :    :   :   :    :   :               : :   :
# : ZS :    : Z :   : Tz :   : -Ry   Rx   D  : : Z :
# 
# 
# Where X,Y,Z are the coordinates in ITRF2008 and XS,YS,ZS are the coordinates in
# the other frames.

# IERS parameters relative to ITRF2008
#
#                 Tx       Ty       Tz       S         Rx        Ry        Rz      date
# ITRF2005       -2.0     -0.9     -4.7      0.94      0.00      0.00      0.00    2000.0
#      rates      0.3      0.0      0.0      0.00      0.00      0.00      0.00
# ITRF2000       -1.9     -1.7    -10.5      1.34      0.00      0.00      0.00    2000.0
#      rates      0.1      0.1     -1.8      0.08      0.00      0.00      0.00
# ITRF97          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF96          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF94          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF93        -24.0      2.4    -38.6      3.41     -1.71     -1.48     -0.30    2000.0
#      rates     -2.8     -0.1     -2.4      0.09     -0.11     -0.19      0.07
# ITRF92         12.8      4.6    -41.2      2.21      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF91         24.8     18.6    -47.2      3.61      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF90         22.8     14.6    -63.2      3.91      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF89         27.8     38.6   -101.2      7.31      0.00      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
# ITRF88         22.8      2.6   -125.2     10.41      0.10      0.00      0.06    2000.0
#      rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
#
# IGS parameters
#
# From 
# Transforming Positions and Velocities between the International Terrestrial Reference Frame of 2000 
# and North American Datum of 1983
# Tomas Soler, M.ASCE, and Richard A. Snay
# JOURNAL OF SURVEYING ENGINEERING (c) ASCE / MAY 2004 / P49
#
# Referenced to 
# Springer, T. A., Kouba, J., and Mireault, Y. 2000. "1999 analysis coor-
# dinator report." 1999 Tech. Rep., International GPS Service for Geo-
# dynamics, Jet Propulsion Laboratory, Pasadena, Calif., 15-55.
#
#  IGS parameters relative to ITRF97
#
#  ITRF96 -2.07 -0.21 9.95 -0.93496 +0.12467 -0.22355 -0.06065 1997.0
#    rates 0.69 -0.10 1.86 -0.19201 +0.01347 -0.01514 +0.00027 
#
# Combined to generate following parameters relative to ITRF96 at epoch 2000.0
#
# Note: IGS parameters use opposite sign convention for rotations to IERS

ITRF_params= (
('ITRF2008', ( -4.8,  -2.09,  17.67, -1.40901, 0.16508, -0.26897, -0.11984),
             ( -0.79,  0.6,    1.34,  0.10201, 0.01347, -0.01514, -0.01973)),
('ITRF2005', ( -6.8,  -2.99,  12.97, -0.46901, 0.16508, -0.26897, -0.11984),
             ( -0.49,  0.6,    1.34,  0.10201, 0.01347, -0.01514, -0.01973)),
('ITRF2000', ( -6.7,  -3.79,   7.17, -0.06901, 0.16508, -0.26897, -0.11984),
             ( -0.69,  0.7,   -0.46,  0.18201, 0.01347, -0.01514, -0.01973)),
('ITRF97',   (  0,     0.51, -15.53,  1.51099, 0.16508, -0.26897, -0.05984),
             ( -0.69,  0.1,   -1.86,  0.19201, 0.01347, -0.01514,  0.00027)),
('ITRF96',   (  0,     0,      0,     0,       0,        0,        0),
             (  0,     0,      0,     0,       0,        0,        0)),
);

ITRF_ref='ITRF96';
refdate=2000.0;

secstorads = math.radians(1.0/3600.0)
scalefactors = (0.001, 0.001, 0.001, 1.0e-9, secstorads*0.001, secstorads*0.001, secstorads*0.001 )

def dateAsYear( dvalue ):
    '''
    Convert a date or datetime to a floating point number of years
    Leaves floating or integer values unchanged.
    '''
    if type(dvalue)==float: 
        return dvalue
    if type(dvalue)==int:
        return float(dvalue)
    if type(dvalue)==datetime:
        year=dvalue.year
        dt=dvalue-datetime(year,1,1)
    elif type(dvalue)==date:
        year=dvalue.year
        dt=dvalue-date(year,1,1)
    else:
        raise RuntimeError(type(dvalue).__name__+' is not a valid type for a date')
    dty=float((date(year+1,1,1)-date(year,1,1)).days)
    return year+(dt.days+dt.seconds/(24.0*60*60))/dty


class BursaWolf14Transformation( object ):

    def __init__( self, rffrom, rfto, params, rates=None, refdate=refdate, source=None ):
        self.rffrom=rffrom
        self.rfto=rfto
        self.params=None if params is None else list(params)
        self.rates=None if rates is None else list(rates)
        self.refdate=refdate
        self.source=source
        self._tf=None

    def __str__(self):
        return ("Transformation from "+self.rffrom+" to "+self.rfto+"\n"+
                (
                " Reference date {0:.1f}\n".format(self.refdate)
                    if self.rates and self.refdate else '') +
                "   Translations {0:.2f}  {1:.2f}  {2:.2f} mm\n".format(*self.params[0:3])+
                (
                "          rates {0:.2f}  {1:.2f}  {2:.2f} mm/yr\n".format(*self.rates[0:3])
                    if self.rates and self.refdate else '') +
                "      Rotations {0:.5f}  {1:.5f}  {2:.5f} mas\n".format(*self.params[4:7])+
                (
                "          rates {0:.5f}  {1:.5f}  {2:.5f} mas/yr\n".format(*self.rates[4:7])
                    if self.rates and self.refdate else '') +
                "          Scale {0:.5f}  ppb\n".format(self.params[3])+
                (
                "          rates {0:.5f} ppb/yr\n".format(self.rates[3])
                    if self.rates and self.refdate else ''))

    def reversed( self ):
        return BursaWolf14Transformation(
            self.rfto,
            self.rffrom,
            [-p for p in self.params],
            None if self.rates is None else [-r for r in self.rates],
            self.refdate,
            self.source
            )

    def atDate( self, date ):
        p=list(self.params)
        refdate=None
        date=dateAsYear(date)
        if self.rates and self.refdate:
            diff = date - self.refdate
            for i,r in enumerate(self.rates):
                p[i]=self.params[i] + r*diff
            refdate = date
        return BursaWolf14Transformation(
            self.rffrom,
            self.rfto,
            p,
            self.rates,
            refdate,
            self.source )

    def add( self, other ):
        if self.rffrom==other.rfto:
            rffrom=other.rffrom
            rfto=self.rfto
        elif self.rfto==other.rffrom:
            rffrom=self.rffrom
            rfto=other.rfto
        else:
            raise RuntimeError("Cannot join incompatible transformations (must have common start/end reference frame")

        refdate=self.refdate if self.refdate is not None else other.refdate
        if refdate and refdate != other.refdate:
            other=other.atDate(refdate)

        return BursaWolf14Transformation(
            rffrom,
            rfto,
            [p1+p2 for p1,p2 in zip(self.params,other.params)],
            (self.rates if other.rates is None 
             else other.rates if self.rates is None 
             else [p1+p2 for p1,p2 in zip(self.rates,other.rates)]),
            refdate,
            self.source if self.source == other.source else None
            )

    def subtract( self, other ):
        return self.add( other.reversed())

    def transFunc( self, date ):
        '''
        Generates a transformation function between the two ITRF
        realisations at a specific date.  The function generated takes
        as parameters a single coordinate or an array of coordinates.

        Uses numpy for efficient processing of arrays of coordinates
        '''
        if self.rates is not None and self.refdate is not None and date is not None:
            diff = dateAsYear(date) - self.refdate
            params=[(p+r*diff)*s for p,r,s in zip(self.params,self.rates,scalefactors)]
        else:
            params=[p*s for p,s in zip(self.params,scalefactors)]
        txyz=np.array([params[0:3]])
        scale=params[3]
        (rx,ry,rz)=params[4:7]
        rxyz=np.transpose(np.array([[scale,-rz,ry],[rz,scale,-rx],[-ry,rx,scale]]))
        def tf( coords ):
            if not isinstance(coords,np.ndarray):
                coords=np.array(coords)
            single=len(coords.shape)==1
            if single:
                coords=coords.reshape((1,coords.size))
            coords=coords+txyz+coords.dot(rxyz)
            if single:
                coords=coords.reshape((coords.size))
            return coords
        return tf

    def transform( self, xyz, date=None ):
        '''
        Transforms a single coordinate [X,Y,Z] or an array of 
        coordinates [[X1,Y1,Z1],[X2,Y2,Z2],..] from the source to the 
        target ITRF.  Transformation is applied at the specified date,
        or at the reference date of the transformation if none is 
        defined.

        Uses numpy for efficient processing of arrays of coordinates
        '''
        if self._tf is None or date != self._tfdate:
            self._tf=self.transFunc(date)
            self._tfdate=date
        return self._tf(xyz)

    def transformLonLat( self, lon, lat=None, hgt=None, date=None ):
        '''
        Transforms a single coordinate [Lon,Lat,Hgt] or an array of 
        coordinates [[lon1,lat1,hgt1],[lon2,lat2,hgt2],..] from the 
        source to the target ITRF.  Transformation is applied at the specified date,
        or at the reference date of the transformation if none is 
        defined.

        Uses numpy for efficient processing of arrays of coordinates
        '''
        from Ellipsoid import GRS80
        xyz=GRS80.xyz(lon,lat,hgt)
        xyz=self.transform(xyz,date=date)
        return GRS80.geodetic(xyz)


def Transformation( to_itrf=ITRF_ref, from_itrf=ITRF_ref ):
    '''
    Determine the transformation from one ITRF realisation to another,
    returns an Transformation object.
    '''
    rffrom=None
    rfto=None
    for p in ITRF_params:
        if p[0] == to_itrf:
            rfto=BursaWolf14Transformation(ITRF_ref,p[0],p[1],p[2],refdate)
            if from_itrf==ITRF_ref:
                return rfto
        if p[0] == from_itrf:
            rffrom=BursaWolf14Transformation(ITRF_ref,p[0],p[1],p[2],refdate).reversed()
            if to_itrf==ITRF_ref:
                return rffrom
    if not rffrom:
        raise RuntimeError(from_itrf+' is not a recognized ITRF')
    if not rfto:
        raise RuntimeError(to_itrf+' is not a recognized ITRF')
    return rffrom.add(rfto)

#itrf2008_nzgd2000=Transformation.transformation(from_itrf='ITRF2008')
#nzgd2000_itrf2008=Transformation.transformation(to_itrf='ITRF2008')

def main():
    import sys
    import argparse
    import re
    parser=argparse.ArgumentParser(description='Convert Cartesian coordinates between ITRF systems')
    parser.add_argument('-f','--from-itrf',default='ITRF2008',help="Source ITRF - default ITRF2008")
    parser.add_argument('-t','--to-itrf',default='ITRF96',help="Target ITRF - default ITRF96")
    parser.add_argument('-d','--date',help="Transformation date (yyyymmdd or yyyy.yyy) - default today")
    parser.add_argument('-l','--list',action='store_true',help="List transformation parameters")
    parser.add_argument('-x','--xyz',nargs=3,type=float,metavar=('X','Y','Z'),help="XYZ coordinates to transform (input/output files ignored)")
    parser.add_argument('-c','--csv',action='store_true',help="File format CSV - default whitespace delimited")
    parser.add_argument('-z','--column-names',metavar=('X_COL','Y_COL','Z_COL'),nargs=3,help="Column names of X,Y,Z fields - default first three columns")
    parser.add_argument('-g','--geodetic',action='store_true',help="Coordinates are lon,lat,hgt")
    parser.add_argument('-v','--verbose',action='store_true',help="More verbose output")
    parser.add_argument('input_file',nargs='?',help="Input file of XYZ coordinates")
    parser.add_argument('output_file',nargs='?',help="Output file of XYZ coordinates")

    args=parser.parse_args()
    input_file=args.input_file
    output_file=args.output_file
    if args.xyz is not None and input_file is not None:
        print "Cannot have xyz and input file arguments"
        sys.exit()
    if not args.list and args.xyz is None and not input_file:
        print "No coordinate input specified - need xyz or input file"
        sys.exit()
    if input_file is not None and output_file is None:
        print "Need the name of an output file"
        sys.exit()

    itrfs=[x[0] for x in ITRF_params]

    from_itrf=args.from_itrf.upper()
    if from_itrf not in itrfs:
        if 'ITRF'+from_itrf in itrfs:
            from_itrf='ITRF'+from_itrf
        else:
            print from_itrf,'is not a valid ITRF'
            print 'Options are:',', '.join(itrfs)
            sys.exit()

    to_itrf=args.to_itrf.upper()
    if to_itrf not in itrfs:
        if 'ITRF'+to_itrf in itrfs:
            to_itrf='ITRF'+to_itrf
        else:
            print to_itrf,'is not a valid ITRF'
            print 'Options are:',', '.join(itrfs)
            sys.exit()

    dt=datetime.now()
    if args.date is not None:
        m = re.match(r'(\d\d\d\d)(\d\d)(\d\d)$',args.date)
        if not m:
            m = re.match(r'(\d\d\d\d)\-(\d\d)\-(\d\d)$',args.date)
        if m:
            dt=datetime(int(m.group(1)),int(m.group(2)),int(m.group(3)))
        elif re.match(r'\d\d\d\d(?:\.\d+)?$',args.date):
            dt=float(args.date)
    year=dateAsYear(dt)

    tfm=Transformation(from_itrf=from_itrf,to_itrf=to_itrf).atDate(year)

    if args.list:
        print tfm
    elif args.verbose:
        print "Transforming from {0} to {1} at {2:.2f}".format(from_itrf,to_itrf,year)

    if args.geodetic:
        transfunc=tfm.transformLonLat
        crdfmt=["{0:.8f}","{0:.8f}","{0:.4f}"]
    else:
        transfunc=tfm.transform
        crdfmt=["{0:.4f}"]*3

    if args.xyz:
        xyzt=transfunc(args.xyz)
        xyzs=[f.format(x) for f,x in zip(crdfmt,xyzt)]
        print "{0} {1} {2}".format(*xyzs)
        sys.exit()

    if args.input_file:
        cols=[0,1,2]
        reqlen=3
        with sys.stdin if input_file=='-' else open(input_file,'r') as fin:
            if args.verbose:
                print "Reading coordinates from",input_file
            import csv
            if args.csv:
                freader=csv.reader(fin)
            else:
                def wsreader(f):
                    for l in f:
                        yield l.split()
                freader=wsreader(fin)
            with sys.stdout if output_file=='-' else open(output_file,'w') as fout:
                if args.verbose:
                    print "Writing coordinates to",output_file
                if args.csv:
                    writerow=csv.writer(fout).writerow
                else:
                    def writerow(row):
                        fout.write('\t'.join(row))
                        fout.write('\n')
                if args.column_names:
                    header=freader.next()
                    writerow(header)
                    header=[x.upper() for x in header]
                    cols=[]
                    for c in args.column_names:
                        c=c.upper()
                        if c not in header:
                            print 'Column',c,'is missing from input file header'
                            sys.exit()
                        cols.append(header.index(c))
                    reqlen=max(cols)
                for row in freader:
                    if len(row) < reqlen:
                        continue
                    xyz=[float(row[i]) for i in cols]
                    xyzt=transfunc(xyz)
                    for f,c,x in zip(crdfmt,cols,xyzt):
                        row[c]=f.format(x)
                    writerow(row)

if __name__ == "__main__":
    main()
