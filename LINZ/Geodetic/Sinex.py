
# Module for reading SINEX files.

import sys
from collections import namedtuple
from datetime import datetime, timedelta
import gzip
import re
import numpy as np

class Reader( object ):

    _scanners={}
    _epochre=re.compile(r'^(\d\d)\:(\d\d\d)\:(\d\d\d\d\d)$')
    
    Site=namedtuple('Site','id code monument description llh')
    Epoch=namedtuple('Epoch','id code soln start end mean')
    Coordinate=namedtuple('Coordinate','id code soln crddate xyz vxyz prmids')
    Solution=namedtuple('Solution','id code soln monument description llh startdate enddate crddate xyz vxyz prmids covariance')

    COVAR_FULL=2
    COVAR_STATION=1
    COVAR_NONE=0

    def __init__( self, filename, **options ):
        '''
        Load a SINEX file. Options are:

            useCodes        a list of codes to load (default is all) 
            readVelocities  read velocities as well as xyz
            covariance      covariance option, one of COVAR_FULL, COVAR_STATION, COVAR_NONE
            readCovariance  read a covariance matrix

        '''
        self._filename=filename
        self._fh=None
        self._lineno=0
        self._version=None
        self._options=options
        self._selectCodes=options.get('selectCodes')
        covarOption=options.get('covariance') or Reader.COVAR_NONE
        if covarOption and covarOption != Reader.COVAR_STATION:
            covarOption=Reader.COVAR_FULL
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
                monuments[monument].append((ptid,ptcode,solnid))
        self._solutions=solutions
        self._monuments=monuments

    def get( self, ptid=None, ptcode=None, solnid=None, monument=None, exceptionIfNone=False ):
        '''
        Get a specific solution.  Can either specify a ptid or a monument.  If this is ambiguous then
        a ptcode and solnid may also be supplied).
        '''
        solutions=[]
        if monument is not None:
            if ptid is not None or ptcode is not None:
                raise RuntimeError('Sinex.Reader.get cannot have monument parameter with ptid or ptcode')
            for solution in self._monuments.get(monument,[]):
                if solnid is not None and solution[2] != solnid:
                    continue
                solutions.append(solution)
        elif ptid is not None:
            for solution in self._solutions.get(ptid,[]):
                if ptcode is not None and solution[1] != ptcode:
                    continue
                if solnid is not None and solution[2] != solnid:
                    continue
                solutions.append(solution)
        else:
            raise RuntimeError('Sinex.Reader.get requires monument or ptid parameter')

        if len(solutions) == 0:
            if exceptionIfNone:
                raise RuntimeError('No Sinex solution found matching request')
            return None
        if len(solutions) > 1:
            raise RuntimeError('Ambiguous Sinex solution requested')

        ptid,ptcode,solnid=solutions[0]
        site=self._sites[ptid,ptcode]
        epoch=self._epochs[ptid,ptcode,solnid]
        coord=self._coords[ptid,ptcode,solnid]
        ncvr,prmids=coord.prmids
        covar=None
        if ncvr < len(self._covariances):
            covar=self._covariances[ncvr]
        else:
            ncvr=None
            prmids=None
        return Reader.Solution(
            ptid,ptcode,solnid,
            site.monument if site is not None else ptid,
            site.description if site is not None else '',
            site.llh if site is not None else None,
            epoch.start if epoch is not None else None,
            epoch.end if epoch is not None else None,
            coord.crddate,
            coord.xyz,
            coord.vxyz,
            prmids,
            covar)

    def _open( self ):
        if self._fh:
            raise RuntimeError('Reopening already open SINEX file handle')
        try:
            if self._filename.endswith('.gz'):
                self._fh=gzip.GzipFile(self._filename)
            else:
                self._fh=open(self._filename)
        except:
            raise RuntimeError('Cannot open SINEX file '+self._filename)
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
                    Reader._scanners[section](self)
                else:
                    self._skipSection(section)
        finally:
            self._close()

    def _readline( self ):
        if self._fh is None:
            raise RuntimeError('Cannot read from unopened SINEX file '+self._filename)
        line=self._fh.readline()
        if line != '':
            self._lineno += 1
        return line.rstrip()

    def _readError( self, message ):
        message=message+' at line '+str(self._lineno)+' of SINEX file '+self._filename
        raise RuntimeError(message)

    def _scanHeader( self ):
        header=self._readline()
        match=re.match(r'^\%\=SNX\s(\d\.\d\d)',header)
        if match is None:
            raise RuntimeError('Invalid SINEX file '+self._filename+' - missing %=SNX in header')

    def _sinexEpoch( self, epochstr ):
        match=self._epochre.match( epochstr )
        if not match:
            self._readError('Invalid epoch string '+epochstr)
        y,d,s=(int(d) for d in match.groups())
        if y==0 and d==0 and s==0:
            return None
        year=y+1900 if y > 50 else y+2000
        return datetime(year,1,1)+timedelta(days=d,seconds=s)

    def _scanSection( self, section, recordre ):
        '''
        Iterates over a section
        '''
        while True: 
            line=self._readline()
            if line == '':
                self._readError(section+' not terminated')
            line=line.rstrip()
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
            if line == '':
                break
            line=line.rstrip()
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
                print line
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

    def _scanSiteId( self ):
        section='SITE/ID'
        recordre=re.compile(r'''^
             \s([\s\w]{4})  # point id
             \s([\s\w]{2})  # point code
             \s([\s\w]{9})  # monument id
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


    def _scanSolutionStatistics( self ):
        section='SOLUTION/STATISTICS'
        recordre=re.compile(r'''^
             \s(.{0,30})     # Statisics item
             \s(.{0,22})     # Statisics value
             \s*$''', re.IGNORECASE | re.VERBOSE )

        stats={}
        for match in self._scanSection(section,recordre):
            stats[match.group(1).strip()]=match.group(2).strip()
        self._stats=stats

    def _scanSolutionEpoch( self ):
        section='SOLUTION/EPOCHS'
        recordre=re.compile(r'''^
             \s([\s\w]{4})  # point id
             \s([\s\w]{2})  # point code
             \s([\s\w\-]{4})  # solution id
             \s(\w)           # to be determined!
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

    def _scanSolutionEstimate( self ):
        section='SOLUTION/ESTIMATE'
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

        usevel=self._options.get('useVelocity')
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
            if (ptid,ptcode,solnid) not in coords:
                prmids=(nextcvr,list(range(nextprm,nextprm+nprm)))
                if covarOption == Reader.COVAR_FULL:
                    nextprm += nprm
                elif covarOption == Reader.COVAR_STATION:
                    nextcvr += 1
                vxyz=np.zeros((3,)) if usevel else None
                coords[ptid,ptcode,solnid]=Reader.Coordinate(ptid,ptcode,solnid,prmtime,
                                                             np.zeros((3,)),vxyz,prmids)
            coord=coords[ptid,ptcode,solnid]
            if prmtime != coord.crddate:
                self._readError('Inconsistent parameter epoch')
            prmlookup[prmid]=(coord.prmids[0],coord.prmids[1][prmno])
            if prmno < 3:
                coord.xyz[prmno]=float(value)
            else:
                coord.vxyz[prmno-3]=float(value)

        self._coords=coords
        self._prmlookup=prmlookup

    def _scanSolutionMatrixEstimate( self ):
        section='SOLUTION/MATRIX_ESTIMATE L COVA'
        recordre=re.compile(r'''^
             \s([\s\d]{5})
             \s([\s\d]{5})
             \s([\s\dE\+\-\.]{21})
             (?:\s([\s\dE\+\-\.]{21}))?
             (?:\s([\s\dE\+\-\.]{21}))?
             \s*$''', re.IGNORECASE | re.VERBOSE )

        if self._covarianceOption == Reader.COVAR_NONE:
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
    }
