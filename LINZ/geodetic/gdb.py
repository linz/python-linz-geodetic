import urllib2
import json
import re
from collections import namedtuple

'''
Module to access information from the LINZ geodetic database.
'''


def _json_object_hook(d): return namedtuple('X',d.keys())(*d.values())

_cache={}
_gdburl='http://www.linz.govt.nz/gdb?mode=js&code={code}'

def get( code, cache=True ):
    '''
    Retrieve information for a geodetic mark. The data is retrieved as an anonymous 
    class (constructed with named tuple) which is built from the JSON returned by the
    geodetic database 'mode=js' option.

    If cache is True then retrieved marks are saved - if the same mark is requested
    again then it is retrieved from the cache.
    '''
    if not re.match(r'^\w{4}$',code):
        raise ValueError(code+' is not a valid geodetic code')
    code=code.upper()
    if cache and code in _cache:
        stn=_cache[code]
    else:
        url=_gdburl.replace('{code}',code)
        try:
            stndata=urllib2.urlopen(url).read()
            stn=json.loads(stndata,object_hook=_json_object_hook)
        except Exception as e:
            raise RuntimeError("Cannot connect to geodetic database: "+e.message)
        if cache:
            _cache[code]=stn
    if stn is None:
        raise ValueError(code+' is not an existing geodetic mark')
    return stn
    
