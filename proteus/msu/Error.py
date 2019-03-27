# Error.py  (burgreen)

# module to provide static error handling

# notes:
# 1. 10oct03 g.burgreen

# from inspect import currentframe, getframeinfo
# frameinfo = getframeinfo(currentframe())
# print( frameinfo.filename, frameinfo.lineno)

import sys

from StringIO import *

#----------------------------------------------
def exit( *args ):
#----------------------------------------------

   s = StringIO()

   if len(args) > 1: 
      str = args[0] % (args[1:])
      print( 'Error.exit: ' + str )
   else: 
      print( 'Error.exit: ' + args[0] )
   
   sys.exit()

#----------------------------------------------
def fatal( *args ):
#----------------------------------------------

   s = StringIO()

   if len(args) > 1:
      #s.write( args[0] % args[1:] )
      #raise RuntimeError, s.getvalue()
      str = args[0] % (args[1:])
      raise RuntimeError, str
   else:
      raise RuntimeError, args[0]

#----------------------------------------------
def warning( *args ):
#----------------------------------------------

   s = StringIO()

   if len(args) > 1:
      #s.write( args[0] % args[1:] )
      #print 'Warning:',s.getvalue()
      str = args[0] % (args[1:])
      print('Error.warning: ' + str )
   else:
      print('Error.warning: ' + args[0] )
