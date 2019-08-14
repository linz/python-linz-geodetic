
# Module for reading SINEX files.

# Imports to support python 3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
from collections import namedtuple
from datetime import datetime, timedelta
import gzip
import math
import re
import numpy as np
from .Ellipsoid import GRS80

try:
    stringtype = basestring   # for Python2
except NameError:
    stringtype = str          # for Python3

COVAR_FULL=2
COVAR_STATION=1
COVAR_NONE=0

EXTRAPOLATE_NONE=0
EXTRAPOLATE_FIRST=1
EXTRAPOLATE_LAST=2
EXTRAPOLATE_AFTER=4
EXTRAPOLATE_ALL=7

_years={}

def _decimalYear( date ):
    global _years
    year=date.year
    if year not in _years:
        y0=datetime(year,1,1)
        y1=datetime(year+1,1,1)
        factor=1.0/(y1-y0).total_seconds()
        _years[year]=(y0,factor)
    return year+(date-_years[year][0]).total_seconds()*_years[year][1]

class NoSolution( RuntimeError ): pass
class AmbiguousSolution( RuntimeError ): pass
class InvalidParameters( RuntimeError ): pass
class CovarianceMissing( RuntimeError ): pass
class SinexFileError( RuntimeError ): pass

class Reader( object ):

    _scanners={}
    _epochre=re.compile(r'^(\d\d)\:(\d\d\d)\:(\d\d\d\d\d)$')
    
    Site=namedtuple('Site','id code monument description llh')
    Epoch=namedtuple('Epoch','id code soln start end mean')
    Coordinate=namedtuple('Coordinate','id code soln crddate xyz vxyz prmids')
    SolutionId=namedtuple('SolutionId','id code soln')
    Solution=namedtuple('Solution','id code soln monument description llh startdate enddate crddate crddate_year xyz vxyz prmids covariance')

    def __init__( self, filename, **options ):
        '''
        Load a SINEX file. Options are:

            selectCodes     a list of codes to load (default is all) 
            velocities      read velocities as well as xyz (boolean, default is True)
            covariance      covariance option, one of COVAR_FULL, COVAR_STATION, COVAR_NONE,
                                default is COVAR_NONE

        '''
        self._filename=filename
        self._fh=None
        self._lineno=0
        self._version=None
        self._options=options
        self._selectCodes=options.get('selectCodes')
        covarOption=options.get('covariance') or COVAR_NONE
        if covarOption and covarOption != COVAR_STATION:
            covarOption=COVAR_FULL
        self._covarianceOption=covarOption
        self._sites={}
        self._stats={}
        self._epochs={}
        self._coords={}
        self._prmlookup={}
        self._covariances=[]
        self._solutionIndex={}
        self._scan()
        self._buildSolutionIndex()

    def solutions( self ):
        '''
        Return a list of solutions as tuples (point id, point code, solution id)
        '''
        return sorted(self._coords.keys())

    def _buildSolutionIndex( self ):
        solutions={}
        monuments={}
        for key in self._coords:
            ptid,ptcode,solnid=key
            if ptid not in solutions:
                solutions[ptid] = []
            solutions[ptid].append((ptid,ptcode,solnid))
            site=self._sites.get((ptid,ptcode))
            if site is not None:
                monument=site.monument
                if monument not in monuments:
                    monuments[monument]=[]
                monuments[monument].append(Reader.SolutionId(ptid,ptcode,solnid))
        self._solutions=solutions
        self._monuments=monuments

    def getSolutionIds( self, ptid=None, ptcode=None, solnid=None, monument=None, solution=None ):
        '''
        Get a list of solution ids matching the point id, point code, and solution id,
        or monument
        or solution tuple

        First parameter may be interpreted as a ptid (stringtype), or ':' separated
        point id, point code, solution id, or solution id tuple (if a tuple), 
        '''

        if isinstance(ptid,tuple):
            solution=ptid
            ptid=None

        if isinstance(solution, stringtype):
            ptid=solution
            solution=None

        if ptid is not None and ':' in ptid:
            solution=((ptid+'::').split(':'))[:3]
            ptid=solution[0]

        if solution is not None:
            ptid=ptid or solution[0]
            ptcode=ptcode or solution[1]
            solnid=solnid or solution[2]

        solutions=[]
        if monument is not None:
            if ptid is not None or ptcode is not None:
                raise InvalidParameters('Sinex.Reader cannot specify monument as well as ptid or ptcode')
            for solution in self._monuments.get(monument,[]):
                if solnid is not None and solution[2] != solnid:
                    continue
                key=(solution[0],solution[1])
                solutions.append(solution)
        elif ptid is not None:
            for solution in self._solutions.get(ptid,[]):
                if ptcode is not None and solution[1] != ptcode:
                    continue
                if solnid is not None and solution[2] != solnid:
                    continue
                solutions.append(solution)
        else:
            raise InvalidParameters('Sinex.Reader requires a solution, monument, ptid parameter')
        return solutions


    def get( self, ptid=None, ptcode=None, solnid=None, monument=None, solution=None, date=None, 
            allSolutions=False, extrapolate=EXTRAPOLATE_NONE ):
        '''
        Get a specific solution.  Can either specify a ptid or a monument.
        If this is ambiguous then a ptcode and solnid may also be supplied). 
        Alternatively can specify a solution such as returned by solutions().

        If allSolutions is true then a list of matching solutions will be returned

        If a date is specified then the solution applicable at that date is returned

        The extrapolate option affects extrapolation beyond the solution dates.  Options are:
            EXTRAPOLATE_NONE    only return solutions matching the date
            EXTRAPOLATE_FIRST   use the first solution for dates before it
            EXTRAPOLATE_LAST    use the last solution for dates after it
            EXTRAPOLATE_AFTER   allow extrapolation each solution forwards till the next 
                                (implies EXTRAPOLATE_LAST)
            EXTRAPOLATE_ALL     use all these options
        Note that extrapolation here refers to the selection of solutions, not the 
        coordinates calculated from them.  Use .xyz to calculate/extrapolate coordinates.
        '''

        solutions=self.getSolutionIds(ptid=ptid,ptcode=ptcode,solnid=solnid,monument=monument,solution=solution)
        if date is not None:

            extrapolateFirst = extrapolate & EXTRAPOLATE_FIRST
            extrapolateAfter = extrapolate & EXTRAPOLATE_AFTER
            extrapolateLast = (extrapolate & EXTRAPOLATE_LAST) or extrapolateAfter

            # Group solution by monument
            msolutions={}
            for s in solutions:
                key=s[:2]
                if key not in msolutions:
                    msolutions[key]=[]
                msolutions[key].append(s)

            # Order by epoch for each monument
            epochs=self._epochs
            for s in msolutions:
                msolutions[s].sort(key=lambda soln: epochs[soln].start if soln in epochs else datetime.min)

            # Select by date
            solutions=[]
            for solnset in msolutions.values():
                # If a date is defined then only return one solution for each monument
                chosen=solnset[0] if extrapolateFirst else None
                for s in solnset:
                    epoch=self._epochs.get(s)
                    if epoch is not None:
                        if not extrapolateAfter:
                            chosen=None
                        if epoch.start > date:
                            break
                        if epoch.start <= date:
                            if epoch.end >= date:
                                chosen=s
                                break
                        if extrapolateLast: # implies extrapolateAfter
                            chosen=s
                if chosen:
                    solutions.append(chosen)

        if len(solutions) == 0:
            raise NoSolution('No Sinex solution found matching request')
            return None
        if len(solutions) > 1 and not allSolutions:
            raise AmbiguousSolution('Ambiguous Sinex solution requested')

        results=self.getSolutions(solutions)
        return results if allSolutions else results[0]

    def getSolutions( self, solutions ):
        '''
        Returns the full solution data for a list of solutions, each supplied
        as a tuple of point id, point code, and solution id (eg as returned by
        solutions, or get).

        If the covariance was loaded (see covariance option of __init__) then
        the solutions will include a covariance matrix.  If the full covariance 
        was loaded each solution will contain the same full covariance matrix.
        Each solution contains prmids array mapping x,y,z and vx, vy, vz into 
        the corresponding elements of the array
        '''
        results=[]
        for solution in solutions:
            ptid,ptcode,solnid=solution[:3]
            site=self._sites.get((ptid,ptcode))
            epoch=self._epochs.get((ptid,ptcode,solnid))
            coord=self._coords[ptid,ptcode,solnid]
            ncvr,prmids=coord.prmids
            covar=None
            if ncvr < len(self._covariances):
                covar=self._covariances[ncvr]
            else:
                ncvr=None
                prmids=None
            results.append(Reader.Solution(
                ptid,ptcode,solnid,
                site.monument if site is not None else ptid,
                site.description if site is not None else '',
                site.llh if site is not None else None,
                epoch.start if epoch is not None else None,
                epoch.end if epoch is not None else None,
                coord.crddate,
                _decimalYear(coord.crddate),
                coord.xyz,
                coord.vxyz,
                prmids,
                covar))
        return results

    def xyz( self, solution=None, date=None, covariance=False, _covarInfo=False,
            extrapolate=EXTRAPOLATE_NONE ):
        '''
        Return the xyz coordinate at an epoch based on the point id and date.  
        Solution can be (ptid,ptcode,solnid) tuple, or "ptid:ptcode:solnid string".
        ptcode can be omitted if ptid is not ambiguous.  solnid can be omitted - 
        will be defined by date.
        '''
        solutions=self.get( solution=solution, date=date,
                           extrapolate=extrapolate, allSolutions=True )

        if not solutions:
            return
        if len(solutions) > 1:
            raise AmbiguousSolution('Ambiguous solution in SINEX file')
        solution=solutions[0]
        if (covariance or _covarInfo) and solution.covariance is None:
            raise CovarianceMissing('Sinex covariances not read')

        if solution.vxyz is not None and date is not None:
            ydiff=_decimalYear(date)-solution.crddate_year
            xyz=np.array(solution.xyz)+ydiff*np.array(solution.vxyz)
            nprm=6
        else:
            ydiff=0.0
            xyz=solution.xyz.copy()
            nprm=3

        if _covarInfo:
            return xyz,ydiff,solution.prmids,nprm,solution.covariance

        if not covariance:
            return xyz

        prmids=solution.prmids
        ptcovar=solution.covariance[np.ix_(prmids,prmids)]
        mult=np.zeros((3,nprm))
        mult[:,0:3]=np.identity(3)
        if nprm > 3:
            mult[:,3:]=np.identity(3)*ydiff
        xyzcovar=mult.dot(ptcovar.dot(mult.T))
        return xyz,xyzcovar

    def dxyz( self, solutionFrom, solutionTo, date=None, enu=False,
             extrapolate=EXTRAPOLATE_NONE, covariance=False):
        '''
        Return the difference in xyz coordinates at an epoch based on the point id and date.  
        Solutions can be supplied as solution tuples or ':' delimited strings of ptid, ptcode, solnid.
        '''

        xyz1=self.xyz(solutionFrom,date,extrapolate=extrapolate,_covarInfo=covariance)
        xyz2=self.xyz(solutionTo,date,extrapolate=extrapolate,_covarInfo=covariance)
        if covariance:
            xyz1,ydiff1,prmids1,nprm1,covar1=xyz1
            xyz2,ydiff2,prmids2,nprm2,covar2=xyz2
            if covar1 is not covar2 or nprm1 != nprm2:
                raise CovarianceMissing('Need full Sinex covariance to calculate vector difference covariance')
            obseq=np.zeros((3,covar1.shape[0]))
            rows=[0,1,2,0,2,1][:nprm1]
            mult=np.ones((nprm1,))
            if nprm1 > 3:
                mult[3:]=ydiff2
            obseq[rows,prmids2[:nprm1]]=mult
            if nprm1 > 3:
                mult[3:]=ydiff1
            obseq[rows,prmids1[:nprm1]] -= mult
            covar=obseq.dot(covar1.dot(obseq.T))

        dxyz=xyz2-xyz1
        if enu:
            xyz=(xyz1+xyz2)/2.0
            llh=GRS80.geodetic(xyz)
            enu_axes=GRS80.enu_axes(llh[0],llh[1])

            dxyz=enu_axes.dot(dxyz)
            if covariance:
                covar=enu_axes.dot(covar.dot(enu_axes.T))
        if covariance:
            return dxyz,covar
        return dxyz

    def _open( self ):
        if self._fh:
            raise SinexFileError('Reopening already open SINEX file handle')
        try:
            if self._filename.endswith('.gz'):
                self._fh=gzip.GzipFile(self._filename)
            else:
                self._fh=open(self._filename)
        except:
            raise SinexFileError('Cannot open SINEX file '+self._filename)
        self._lineno=0
        self._scanHeader()

    def _close( self ):
        self._fh.close()
        self._fh=None
        self._lineno=0

    def _scan( self ):
        self._open()
        sections=[]
        try:
            while True:
                section=self._nextSection()
                if section is None:
                    break
                sections.append(section)
                if section in Reader._scanners:
                    Reader._scanners[section](self,section)
                    if section=='SOLUTION/ESTIMATE' and self._covarianceOption==COVAR_NONE:
                        break
                else:
                    self._skipSection(section)
        finally:
            self._close()

    def _readline( self ):
        if self._fh is None:
            raise SinexFileError('Cannot read from unopened SINEX file '+self._filename)
        line=self._fh.readline()
        if line == '': return None
        self._lineno += 1
        return line.rstrip()

    def _readError( self, message ):
        message=message+' at line '+str(self._lineno)+' of SINEX file '+self._filename
        raise SinexFileError(message)

    def _scanHeader( self ):
        header=self._readline()
        match=re.match(r'^\%\=SNX\s(\d\.\d\d)',header)
        if match is None:
            raise SinexFileError('Invalid SINEX file '+self._filename+' - missing %=SNX in header')

    def _sinexEpoch( self, epochstr ):
        match=self._epochre.match( epochstr )
        if not match:
            self._readError('Invalid epoch string '+epochstr)
        y,d,s=(int(d) for d in match.groups())
        if y==0 and d==0 and s==0:
            return None
        year=y+1900 if y > 50 else y+2000
        return datetime(year,1,1)+timedelta(days=d-1,seconds=s)

    def _scanSection( self, section, recordre ):
        '''
        Iterates over a section
        '''
        while True: 
            line=self._readline()
            if line is None:
                self._readError(section+' not terminated')
            if line == '':
                continue
            cmdchar=line[0]
            if cmdchar=='*':
                continue
            if cmdchar=='-':
                if line[1:].strip() != section:
                    self._readError(section+' terminated incorrectly')
                break
            if cmdchar != ' ':
                self._readError(section+' invalid line')
            if recordre is not None:
                match=recordre.match(line)
                if match is None:
                    self._readError(section+' badly formatted line: '+line)
                yield match

    def _skipSection( self, section ):
        for match in self._scanSection( section, None ):
            pass

    def _nextSection( self ):
        section=None
        while True:
            line=self._readline()
            if line is None:
                break
            if line == '%ENDSNX':
                break
            if line=='':
                continue
            cmdchar=line[0]
            if cmdchar == '*':
                continue
            if cmdchar == ' ':
                continue
            if cmdchar != '+':
                self._readError('Unexpected character '+cmdchar+' looking for section')
            section=line[1:]
            break
        return section

    def _wantId( self, ptid ):
        return self._selectCodes is None or ptid in self._selectCodes

    def _latlon( self, llstring ):
        negative='-' in llstring
        llstring.replace('-',' ')
        d,m,s=(float(x) for x in llstring.split())
        angle=d+m/60.0+s/3600.0
        if negative:
            angle=-angle
        return angle

    def _scanSiteId( self, section ):
        recordre=re.compile(r'''^
             \s([\s\w]{4})  # point id
             \s([\s\w]{2})  # point code
             \s(.{9})  # monument id
             \s[\w\s]       # Observation techniquesl
             \s(.{22})      # description
             \s([\s\-\d]{3}[\-\s\d]{3}[\-\s\d\.]{5}) # longitude DMS
             \s([\s\-\d]{3}[\-\s\d]{3}[\-\s\d\.]{5}) # latitude DMS
             \s([\s\d\.\-]{7})                       # height
             \s*$''', re.IGNORECASE | re.VERBOSE )

        sites={}
        for match in self._scanSection(section,recordre):
            ptid,ptcode,ptname,description,lon,lat,hgt=(x.strip() for x in match.groups())
            if self._wantId( ptid ):
                llh=(self._latlon(lon),self._latlon(lat),float(hgt))
                sites[ptid,ptcode]=Reader.Site(ptid,ptcode,ptname,description,llh)
        self._sites=sites


    def _scanSolutionStatistics( self, section ):
        recordre=re.compile(r'''^
             \s(.{0,30})     # Statisics item
             \s(.{0,22})     # Statisics value
             \s*$''', re.IGNORECASE | re.VERBOSE )

        stats={}
        for match in self._scanSection(section,recordre):
            stats[match.group(1).strip()]=match.group(2).strip()
        self._stats=stats

    def _scanSolutionEpoch( self, section ):
        recordre=re.compile(r'''^
             \s([\s\w]{4})  # point id
             \s([\s\w]{2})  # point code
             \s([\s\w\-]{4})  # solution id
             \s([\w|\s])           # could be blank
             \s(\d\d\:\d\d\d\:\d\d\d\d\d) # start epoch
             \s(\d\d\:\d\d\d\:\d\d\d\d\d) # end epoch
             \s(\d\d\:\d\d\d\:\d\d\d\d\d) # mean epoch
             \s*$''', re.IGNORECASE | re.VERBOSE )

        epochs={}
        for match in self._scanSection(section,recordre):
            ptid,ptcode,solnid,skip,start,end,mean=(x.strip() for x in match.groups())
            starttime=self._sinexEpoch(start)
            endtime=self._sinexEpoch(end)
            meantime=self._sinexEpoch(mean)
            epochs[ptid,ptcode,solnid]=Reader.Epoch(ptid,ptcode,solnid,starttime,endtime,meantime)
        self._epochs=epochs

    def _scanSolutionEstimate( self, section ):
        recordre=re.compile(r'''^
             \s([\s\d]{5})  # param id
             \s([\s\w]{6})  # param type
             \s([\s\w]{4}|[\s\-]{4})  # point id
             \s([\s\w]{2}|[\s\-]{2})  # point code
             \s([\s\w]{4}|[\s\-]{4})  # solution id
             \s(\d\d\:\d\d\d\:\d\d\d\d\d) #parameter epoch
             \s([\s\w\/]{4})  # param units
             \s([\s\w]{1})  # param constraints
             \s([\s\dE\+\-\.]{21})  # param value
             \s([\s\dE\+\-\.]{11})  # param stddev
             \s*$''', re.IGNORECASE | re.VERBOSE )

        coords={}
        prmlookup={}

        usevel=self._options.get('velocities',True)
        useprms=('STAX','STAY','STAZ','VELX','VELY','VELZ') if usevel else ('STAX','STAY','STAZ')
        nprm=6 if usevel else 3
        covarOption=self._covarianceOption
        nextprm=0
        nextcvr=0

        for match in self._scanSection(section,recordre):
            prmid,prmtype,ptid,ptcode,solnid,epoch,units,constraint,value,stddev=(
                x.strip() for x in match.groups())
            if not self._wantId( ptid ):
                continue
            if prmtype not in useprms:
                continue
            prmno=useprms.index(prmtype)
            prmtime=self._sinexEpoch(epoch)
            prmid=int(prmid)
            solnid=Reader.SolutionId(ptid,ptcode,solnid)
            if solnid not in coords:
                prmids=(nextcvr,(range(nextprm,nextprm+nprm)))
                if covarOption == COVAR_FULL:
                    nextprm += nprm
                elif covarOption == COVAR_STATION:
                    nextcvr += 1
                vxyz=np.zeros((3,)) if usevel else None
                coord=Reader.Coordinate(ptid,ptcode,solnid,prmtime,np.zeros((3,)),vxyz,prmids)
                coords[solnid]=coord
            coord=coords[solnid]
            if prmtime != coord.crddate:
                self._readError('Inconsistent parameter epoch')
            prmlookup[prmid]=(coord.prmids[0],coord.prmids[1][prmno])
            if prmno < 3:
                coord.xyz[prmno]=float(value)
            else:
                coord.vxyz[prmno-3]=float(value)

        self._coords=coords
        self._prmlookup=prmlookup

    def _scanSolutionMatrixEstimate( self, section ):
        recordre=re.compile(r'''^
             \s([\s\d]{5})
             \s([\s\d]{5})
             \s([\s\dE\+\-\.]{21})
             (?:\s([\s\dE\+\-\.]{21}))?
             (?:\s([\s\dE\+\-\.]{21}))?
             \s*$''', re.IGNORECASE | re.VERBOSE )

        if self._covarianceOption == COVAR_NONE:
            self._skipSection( section )
            return

        ncvr=0
        nprm=0
        for cvr,prm in self._prmlookup.values():
            if cvr > ncvr:
                ncvr=cvr
            if prm > nprm:
                nprm=prm
        ncvr+=1
        nprm+=1
        covariances=[]
        for i in range(ncvr):
            covariances.append(np.zeros((nprm,nprm)))
        prmlookup=self._prmlookup


        for match in self._scanSection(section,recordre):
            prmid1=int(match.group(1).strip())
            if prmid1 not in prmlookup:
                continue
            prms1=prmlookup[prmid1]
            ncvr=prms1[0]
            nprm1=prms1[1]

            prmid2=int(match.group(2).strip())
            for i in range(3):
                cvrval=match.group(3+i)
                if cvrval is None:
                    break
                cvrval=cvrval.strip()
                if cvrval == '':
                    break
                if prmid2 in prmlookup:
                    prms2=prmlookup[prmid2]
                    if prms2[0]==ncvr:
                        nprm2=prms2[1]
                        cvrval=float(cvrval)
                        covariances[ncvr][nprm1,nprm2]=cvrval
                        if nprm1 != nprm2:
                            covariances[ncvr][nprm2,nprm1]=cvrval
                prmid2 += 1

        self._covariances=covariances


Reader._scanners={
        'SITE/ID': Reader._scanSiteId,
        'SOLUTION/STATISTICS': Reader._scanSolutionStatistics,
        'SOLUTION/EPOCHS': Reader._scanSolutionEpoch,
        'SOLUTION/ESTIMATE': Reader._scanSolutionEstimate,
        'SOLUTION/MATRIX_ESTIMATE L COVA': Reader._scanSolutionMatrixEstimate,
        'SOLUTION/MATRIX_ESTIMATE U COVA': Reader._scanSolutionMatrixEstimate,
        }


class Writer( object ):
    '''
    The Sinex.Writer class is most simply used as a context manager.  The usage is
    with Sinex.Writer(filename) as snx:
        snx.addFileInfo(...)
        snx.setObsDateRange(startdate,enddate)
        snx.setSolutionStatistics(varianceFactor,sumSquaredResiduals,nObs,nUnknowns )
        for m in marks:
            snx.addMark(m.code, ... )
            snx.addSolution(m,xyz,xyzprms)
        snx.setCovariance(covar)
    Alternatively can be used as an ordinary object, in which case a final call
    to snx.write() is required to create the file.
    '''
    # For the future - make consistent with Reader, make a Sinex object that can be 
    # processed after reading or constructing...

    Mark=namedtuple('Mark','mark id code monument description solutions')
    Solution=namedtuple('Solution','solnid xyz params vxyz vparams startdate enddate crddate')

    def __init__( self, filename, **info ):
        self._filename=filename
        self._fileRefInfo={}
        self._comments=[]
        self._written=False
        self._agency='---'
        self._version='2.10'
        self._constraint='2'
        self._stats=None
        self._startdate=None
        self._enddate=None
        self._marks={}
        self._covariance=None
        self._varfactor=1.0
        self.addFileInfo( **info )

    def __enter__( self ):
        return self

    def __exit__(self,exctype,value,traceback):
        if exctype is None:
            if not self._written:
                self.write()

    def addFileInfo( self, **info ):
        '''
        Add file header information, any of:
            description
            output
            contact
            software
            hardware
            input
            agency
            comment
            constraint
        '''
        for item, value in info.iteritems():
            item=item.upper()
            if item in ('DESCRIPTION','OUTPUT','CONTACT','SOFTWARE','HARDWARE','INPUT'):
                self._fileRefInfo[item]=value
            elif item == 'COMMENT':
                self._comments.extend(value.split('\n'))
            elif item == 'AGENCY':
                self._agency=value.upper()[:3]
            elif item == 'CONSTRAINT':
                if value not in ('0','1','2'):
                    raise SinexFileError('Invalid SINEX constraint code '+value+' specified')
                self._constraint=value

    def setSolutionStatistics( self, varianceFactor, sumSquaredResiduals, nObs, nUnknowns ):
        '''
        Set the solution statistics values
           variance factor
           sum of squared residuals
           number of observations
           number of unknowns
        '''
        self.stats=[varianceFactor,sumSquaredResiduals,nObs, nUnknowns ]

    def setObsDateRange( self, startdate, enddate ):
        '''
        Set the date range of the observations
        '''
        self._startdate=startdate
        self._enddate=enddate

    def addMark( self, mark, id=None, code=None, monument=None, description=None ):
        '''
        Add a mark to the solution.  Can specify the id, code, monument name, and
        description used to identify the mark
        '''
        id=(id or mark).upper()[:4]
        code=(code or 'A')[:2]
        monument=(monument or mark).upper()[:9]
        description=description or ''
        self._marks[mark]=Writer.Mark(mark,id,code,monument,description,[])

    def addSolution( self, mark, xyz, params, vxyz=None, vparams=None, startdate=None, enddate=None, crddate=None):
        '''
        Add a coordinate solution to the sinex file.  Specifies:
            mark (which must have been added with addMark), 
            xyz coordinate components [x,y,z], 
            params efining the row numbers of each ordinate in the covariance matrix (px,py,pz),
            vxyz (optional) velocity components
            vparams (optional) parameter numbers of velocity components
            startdate (optional) start date of applicability of the solution
            enddate (optional) end date of applicability of the solution
        '''
        if mark not in self._marks:
            raise SinexFileError("Cannot add solution for "+mark+" to SINEX file: mark not defined")
        if len(xyz) != 3 or len(params) != 3:
            raise SinexFileError("Invalid coordinate xyz or params specified for mark "+mark)
        params=list(params)
        if vxyz is not None:
            if len(vxyz) != 3 or len(vparams) != 3:
                raise SinexFileError("Invalid coordinate xyz or params specified for mark "+mark)
            vparams=list(vparams)
        else:
            vparams=None
        
        solnid="{0:04d}".format(len(self._marks[mark].solutions)+1)
        self._marks[mark].solutions.append(
            Writer.Solution(solnid,xyz,params,vxyz,vparams,startdate,enddate,crddate))

    def setCovariance( self, covariance, varianceFactor=1.0 ):
        '''
        Define the covariance matrix for the ordinates. (numpy nxn array)
        '''
        self._covariance=covariance
        self._varfactor=varianceFactor

    def _sinexEpoch( self, date=None  ):
        if date is None:
            date=datetime.now()
        diff=date-datetime(date.year,1,1)
        year=date.year-1900
        if year > 99:
            year -= 100
        return "{0:02d}:{1:03d}:{2:05d}".format(year,diff.days+1,diff.seconds)

    def _sinexDms( self, angle ):
        sign=-1 if angle < 0 else 1
        angle=abs(angle)
        angle += 1.0/(3600*10) # For rounding seconds to 1dp
        deg=int(angle)
        angle = (angle - deg)*60
        min=int(angle)
        angle=(angle-min)*60
        angle -= 0.1 # Finish rounding fix
        if angle < 0:
            angle=0.0
        if deg > 0:
            deg *= sign
        elif min > 0:
            min *= sign
        else:
            angle *= sign
        return "{0:3d}{1:3d}{2:5.1f}".format(deg,min,angle)


    def write( self ):
        '''
        Write the file - can only be called once.  If Writer is used as a context manager
        this is called when the context closes.
        '''
        if self._written:
            raise SinexFileError("Cannot write already written SINEX file")
        if len(self._covariance) is None:
            raise SinexFileError('Cannot write SINEX file: no covariance matrix defined')
        # After this point supplied data may be modified!
        self._written=True
        marks=[]
        solutions=[]
        covarprms=[-1]
        marklist={}
        prmid=0
        for m in sorted(self._marks.values(),key=lambda m: (m.id,m.code,m.monument)):
            if len(m.solutions) == 0:
                continue
            key=(m.id,m.code,m.monument)
            if key in marklist:
                raise SinexFileError('Duplicate mark ({0} {1} {2} in SINEX file'.format(*key))
            marklist[key]=1
            marks.append(m)
            for s in m.solutions:
                for i,p in enumerate(s.params):
                    covarprms.append(p)
                    prmid += 1
                    s.params[i]=prmid
                if s.vparams is not None:
                    for i,p in enumerate(s.params):
                        covarprms.append(p)
                        prmid += 1
                        s.params[i]=prmid

        if prmid == 0:
            raise SinexFileError('Cannot write SINEX file: no solution coordinates defined')
        
        startdate=self._startdate or datetime.now()
        enddate=self._startdate or datetime.now()
        varf=self._varfactor

        with open(self._filename,'w') as sf:
            sf.write("%=SNX {0:4.4s} {1:3.3s} {2} {1:3.3s} {3} {4} C {5:5d} {6} S\n"
                     .format(self._version,self._agency,
                             self._sinexEpoch(),
                             self._sinexEpoch(startdate),
                             self._sinexEpoch(enddate),
                             prmid,
                             self._constraint))

            sf.write("+FILE/REFERENCE\n")
            for item in ('DESCRIPTION','OUTPUT','CONTACT','SOFTWARE','HARDWARE','INPUT'):
                sf.write(" {0:18.18s} {1:.60s}\n".format(
                    item,
                    self._fileRefInfo.get(item,'') or 'N/A'))
            sf.write("-FILE/REFERENCE\n")

            if len(self._comments) > 0:
                sf.write("+FILE/COMMENT\n")
                for c in self.comments:
                    sf.write(" {0:.79s}\n".format(c))
                sf.write("-FILE/COMMENT\n")

            sf.write("+SITE/ID\n")
            for m in marks:
                xyz=m.solutions[0].xyz
                lon,lat,hgt=GRS80.geodetic(xyz)
                if lon < 0:
                    lon += 360
                sf.write(" {0:4.4s} {1:2.2s} {2:9.9s} C {3:22.22s} {4} {5} {6:7.1f}\n"
                         .format(m.id,m.code,m.monument,m.description,
                                 self._sinexDms(lon), self._sinexDms(lat),hgt))

            sf.write("-SITE/ID\n")

            # sf.write("+SITE/DATA\n")
            # sf.write("-SITE/DATA\n")

            sf.write("+SOLUTION/EPOCHS\n")
            for m in marks:
                for s in m.solutions:
                    startdate=s.startdate or startdate
                    enddate=s.enddate or enddate
                    diff=enddate-startdate
                    middate=startdate+timedelta(seconds=diff.total_seconds()/2)
                    sf.write(" {0:4.4s} {1:2.2s} {2:4.4s} C {3} {4} {5}\n".format(
                             m.id,m.code,s.solnid,
                             self._sinexEpoch(startdate),
                             self._sinexEpoch(enddate),
                             self._sinexEpoch(middate),
                            ))
            sf.write("-SOLUTION/EPOCHS\n")

            if self.stats is not None:
                sf.write("+SOLUTION/STATISTICS\n")
                sf.write(" {0:30.30s} {1:22.15e}\n".format(
                    'VARIANCE FACTOR',self.stats[0]))
                sf.write(" {0:30.30s} {1:22.15e}\n".format(
                    'SUM OR SQUARED RESIDUALS',self.stats[1]))
                sf.write(" {0:30.30s} {1:22d}\n".format(
                    'NUMBER OF OBSERVATIONS',self.stats[2]))
                sf.write(" {0:30.30s} {1:22d}\n".format(
                    'NUMBER OF UNKNOWNS',self.stats[3]))
                sf.write("-SOLUTION/STATISTICS\n")

            sf.write("+SOLUTION/ESTIMATE\n")
            crds=('X','Y','Z')
            for m in marks:
                for s in m.solutions:
                    startdate=s.startdate or startdate
                    enddate=s.enddate or enddate
                    diff=enddate-startdate
                    crddate=s.crddate or startdate+timedelta(seconds=diff.total_seconds()/2)
                    for crd,value,prmno in zip(crds,s.xyz,s.params):
                        cvrprm=covarprms[prmno]
                        stderr=math.sqrt(varf*self._covariance[cvrprm,cvrprm])
                        sf.write(" {0:5d} STA{1:1.1s}   {2:4.4s} {3:2.2s} {4:4.4s} {5} {6:4.4s} {7:1.1s} {8:21.14e} {9:11.5e}\n"
                                 .format(prmno,crd,m.id,m.code,s.solnid,
                                         self._sinexEpoch(crddate),'m','2',
                                         value,stderr))

            sf.write("-SOLUTION/ESTIMATE\n")

            sf.write("+SOLUTION/MATRIX_ESTIMATE L COVA\n")
            for i in range(1,prmid+1):
                ci=covarprms[i]
                for j in range(1,i+1):
                    if j % 3 == 1:
                        if j > 1:
                            sf.write("\n")
                        sf.write(" {0:5d} {1:5d}".format(i,j))
                    cj=covarprms[j]
                    sf.write(" {0:21.14e}".format(varf*self._covariance[ci,cj]))
                sf.write("\n")
            sf.write("-SOLUTION/MATRIX_ESTIMATE L COVA\n")

            sf.write("%ENDSNX\n")
