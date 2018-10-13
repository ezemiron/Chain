Acquisition and analysis of SIM data
=============================================

Where are the figures?
----------------------

On a fresh clone of the repository, the `figures` repository will be empty.
This is by design.  We should be able to generate the figures in a
reproducible manner from our raw data and scripts.  As part of the build
process, the raw data is re-analyzed and fresh figures are generated.


Where is the raw data?
----------------------

The raw data is far too large to keep under version control (each SIM result
is around 150MB after conversion to 16bit).  At the moment it is assumed that
the files exist in the `data` directory but we will come up with something
better later.  Maybe using our omero server to serve this files?

