import os

if 'HOSTNAME' in os.environ and (os.environ['HOSTNAME'].startswith('garnet')
                                 or
                                 os.environ['HOSTNAME'].startswith('copper')):
    from garnet import *
elif 'HOSTNAME' in os.environ and os.environ['HOSTNAME'].startswith('viutill'):
    from viutill import *
else:
    from default import *