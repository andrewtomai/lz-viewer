##setup.py

from setuptools import setup
setup(
  name = 'lz_assoc_viewer',
  packages = ['lz_assoc_viewer'], # this must be the same as the name above
  version = '1.0',
  description = 'Locuszoom genome wide association viewer plots',
  author = 'Andrew Tomai',
  author_email = 'andrewtomai23@gmail.com',
  license = "MIT",
  url = 'https://github.com/atomai/lz_viewer', # use the URL to the github repo
  download_url = 'https://github.com/atomai/lz_viewer/tarball/1.0', # I'll explain this in a second
  keywords = ['LocusZoom ', 'plot ', 'association ', 'genome '], # arbitrary keywords
  package_data={
	  '': ['static/*'], 
	  '': ['templates/*'],
	  },
  classifiers = [],
)