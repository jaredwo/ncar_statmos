#############################################################################################
NCAR STATMOS/SAMSI Project: Representing Extremes in High-Resolution Gridded Climate Products
#############################################################################################

This project is part of the `NCAR STATMOS/SAMSI workshop <https://sites.google.com/a/uchicago.edu/ncar17>`_ held at NCAR from July
17-21, 2017. The objective of the project is to investigate spatial interpolation
methods for climate data that eliminate or minimize smoothing, so that extreme values
are better represented. Original project proposal can be found `here <https://drive.google.com/open?id=1N8NWwQhLh1LSvf6RuMl4B-VTulDqPBxo0kxo3KYRUAE>`_.

Project participants
====================

* Liz Drenkard, Scripps Institution of Oceanography, liz.drenkard@gmail.com
* Hossein Moradi Rekabdarkolaee, Virginia Commonwealth University, moradirekabh@vcu.edu
* Jared Oyler, Penn State University, jared.oyler@psu.edu

Data
====================
The main analysis data for the project are 1948-2015 daily precipitation observations from
the `GHCN-D archive <https://www.ncdc.noaa.gov/ghcn-daily-description>`_. Spatial
domains of analysis include the U.S. Mid-Atlantic and Colorado. Data can be found
on the `Google Drive folder <https://drive.google.com/open?id=0B9TBb2nhU2d9LWYzODRMT2phbE0>`_ for the project.

Data are provided in `netcdf format <https://www.unidata.ucar.edu/software/netcdf/>`_.
Observations are stored as 2D space-wide matrices where each row is a day in time and
each column is a station. Missing values are represented as nan. Examples for reading
the netcdf files are provided for both R and Python.

Background reading
====================

* `High resolution gridded climate products <https://paperpile.com/shared/RpFjxd>`_
* `Global vs. local accuracy of interpolation methods <https://paperpile.com/shared/ts4lfT>`_

