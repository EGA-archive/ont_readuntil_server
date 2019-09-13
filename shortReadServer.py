# shortReadServer.py

# Index files are created with BWA. Specify the reference file name in 'index'
# The Mask file is 1 byte per position, divided into upper and lower 4 bits
# for forward and reverse strands.
# Specify the the name of the map file in 'mapfile'

# Usage: curl 127.0.0.1:8000/map?TTTACG
# Specify the sequence as query flag "?TTTACG" with URL "map" on port 8000
# The server will return '-1' to indicate stop or '1' to indicate continue
# '0' to indicate wait has not been assigned yet
# (The algorithm to properly assign this vaue is still TODO)

import os
import falcon
import numpy as np
from bwapy import BwaAligner

class SRS:
    index = '/home/asenf/devel/data/aae_ref_AaegL5.0_chr1.fa.gz'
    mapfile = 'data.dat'

    def __init__(self):  
        print('Loading Index: ' + self.index)
        #options = '-x ont2d -A 1 -B 0'
        #self.aligner = BwaAligner(self.index, options=options)
        self.aligner = BwaAligner(self.index)
        print('Index Loaded')
        print('Accessing Mapping File: ' + self.mapfile)
        size = os.path.getsize(self.mapfile)
        print('Mask file is size: ' + str(size))
        self.fp = np.memmap(self.mapfile, dtype='int8', mode='readonly', shape=(1,size))
        print('Map File Accessed')

    def on_get(self, req, resp):
        response_code = -1
        hits = 0
        hits_on_target = 0
        hits_off_target = 0
        alignments = self.aligner.align_seq(req.query_string)
        for aln in alignments:
            val = self.fp[0][aln.pos]
            val_low = np.bitwise_and(val, 0x0f)
            val_up = np.bitwise_and(val >> 4, 0x0f)
            hits += 1
            if val_up > 0 or val_low > 0:
                response_code = 1
                hits_on_target += 1
            else:
                hits_off_target += 1

        resp.media = {'Answer':response_code,
                      'Hits':hits,
                      'HitsOnTarget':hits_on_target,
                      'HitsOffTarget':hits_off_target}

api = falcon.API()
api.add_route('/map', SRS())
