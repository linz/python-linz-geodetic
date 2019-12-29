import sys
import os
import os.path
import re
import unittest
import numpy as np
from numpy import array

# Set output options for dumping numpy data
np.set_printoptions(precision=10,suppress=True)

# Floating point test tolerance
defaultDelta=1.0e-8

# Set True to write the results from all tests to stdout
dumpResults=False
dumpFileReset=False
try:
    stringtype = str
except NameError:
    stringtype = str

class TestCase( unittest.TestCase ):
    '''
    Subclass of unittest.TestCase to support generation and use of a file
    of test results.  Basic structure of a test module using this is:

        from LINZ import fileunittest

        class MyTestCase( fileunittest.TestCase ):

            def test_001( self ):
                result=calc_result()
                self.check('test 001 result',result)

            def test002( self ):
                if notgood():
                    self.reportFail('Test 2 failed')

        if __name__ == "__main__":
            fileunittest.main()

    Run with the --dump option to create a file with the same name as the 
    test module but extension .results.new  Run without --dump to test against
    file .results.
    '''

    resultsFile=None
    dumpFile=None
    testResults={}

    @classmethod
    def setUpClass( cls ):
        import inspect
        clsfile=inspect.getfile(cls)
        TestCase.resultsFile=os.path.splitext(clsfile)[0]+'.results'
        TestCase.dumpFile=os.path.splitext(clsfile)[0]+'.results.new'
        data=''
        try:
            with open(TestCase.resultsFile) as rf:
                data=rf.read()
        except:
            pass
        resultRe=re.compile(r'^\>\>\>\>\s+(\w+)\s*(.*?)(?=^\>\>\>\>)',re.MULTILINE | re.DOTALL)
        for test,result in resultRe.findall(data):
            TestCase.testResults[test]=result.strip()

    def setUp( self ):
        global dumpResults
        global dumpFileReset
        self.dumpfh=None
        if dumpResults: 
            if not dumpFileReset:
                if os.path.exists(self.dumpFile):
                    os.unlink(self.dumpFile)
                dumpFileReset=True
            self.dumpfh=open(self.dumpFile,'a')
            self.dumpfh.write(("="*70)+"\n");

    def tearDown( self ):
        if self.dumpfh:
            self.dumpfh.write(">>>>\n")
            self.dumpfh.close()

    def check( self,testname,output,message=None,delta=None):
        '''
        Call to check that output matches the value for testname in the .results
        file.  message is the error message if it fails.  delta is a numerical
        tolerance for floating point tests.  delta is used in the unittest.assertAlmostEqual
        call.
        '''
        global dumpResults
        testcode=testname.lower()
        testcode=re.sub(r'\s+','_',testcode)
        testcode=re.sub(r'\W','_',testcode)
        if dumpResults:
            self.dumpfh.write(">>>> "+testcode+" ")
            self.dumpfh.write(output if isinstance(output,stringtype) else repr(output))
            self.dumpfh.write("\n")
        else:
            message=message or testname+' incorrect'
            expected=self.testResults.get(testcode)
            if not isinstance(output,str):
                try:
                    expected=eval(expected)
                except:
                    output=repr(output)
            self.checkEqual(output,expected,message,delta)

    def checkRun( self, testname, function, message=None, delta=None ):
        # Run a function with no parameters and either check the output
        # or check the error generated.
        try:
            self.check(testname,function(),message,delta)
        except Exception as ex:
            self.check(testname,ex,message,delta)

    def checkEqual( self, output, expected, message, delta ):
        # Equality check with tolerance on floating point numbers and handling
        # of numpy arrays
        global defaultDelta
        delta=delta or defaultDelta
        if isinstance(output,np.ndarray):
            error=np.max(np.abs(output-expected))
            if error > delta:
                self.fail(message+" max diff {0:.3e}".format(error))
        elif isinstance(output,float):
            message=message+" ({0} != {1})".format(output,expected)
            self.assertAlmostEqual(output,expected,msg=message,delta=delta)
        elif isinstance(output,str):
            output=output.strip()
            if isinstance(expected,str):
                expected=expected.strip()
            if output != expected:
                self.fail(message+"  ({0} != {1})".format(output,expected))
        elif isinstance(output,dict):
            if not message.endswith(']'):
                message=message+' '
            for k in list(expected.keys()):
                self.checkEqual(output[k],expected[k],message+'[{0}]'.format(k),delta)
        elif hasattr(output,'__getitem__'):
            if not message.endswith(']'):
                message=message+' '
            for i,e in enumerate(expected):
                o=output[i]
                self.checkEqual(o,e,message+'[{0}]'.format(i),delta)
        else:
            self.assertEqual(output,expected,msg=message)

    def reportFail( self, message ):
        '''
        Function to directly report a failed test
        '''
        global dumpResults
        if not dumpResults:
            self.fail(message)

def main():
    '''
    Main function called from sublclasses to run the tests directly
    '''
    global dumpResults
    global dumpFileReset
    if '--dump' in sys.argv:
        sys.argv.remove('--dump')
        dumpResults=True
        dumpFileReset=False
        print("**** Dumping output - not running tests!")
    unittest.main()
