# ont_readuntil_server
A Python-based Short Loop / Short Read Masked Match Server

This project is implemented as a gunicorn/falcon based REST server. 

Dependencies: gunicorn, falcon, pybwa, numpy

To start type: gunicorn --workers=1 --threads=4 shortReadServer:api

The parameter '--workers' starts independent threads of the server code, this allows effective scaling of server performance, but it is memory-intensive, because it duplicates the in-memnory footprint of the BWA index file.
The parameter '--threads' increases the parallel threads within the Python environment; but this has little positive effect on performance, due to Python thread handling issues.

The server implements two endpoints:

../map?{seq} [where {seq} is a short sequence string]
../index

A call to /index re-loads the index file used by bwa into memory. This functionality is necessary because it is anticipated that the reference sequence will change in the course of a sequencing run, and to correctly decide which sequences to reject via Read Until, it is necessary to update the mask + reference files according to what is observed.

A call to /map along with a gene sequence performs these steps: (1) map the seqence to the reference, using pybwa. (2) look at the mask value at that position, and return an answer corresponding to "continue sequencing" and "stop".

The full return is: 
                     {'Answer':response_code,
                      'Hits':hits,
                      'HitsOnTarget':hits_on_target,
                      'HitsOffTarget':hits_off_target}

'Hits' total number of matches found for a sequence. 'OnTarget' Count of 'continue' answers according to tha mask; 'OffTarget' count of 'stop' answers. The response code is "continue" is at least one match produces this answer.

The location (paths) of the index and the map files are hard-coded in the Python script, and must be adjusted to the local environment. 
