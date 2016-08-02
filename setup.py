from setuptools import setup

setup(name='geo_frac_analysis',
      version='0.1.4',
      description='Python module for the spatial analysis of geological fractures.',
      long_description='Python module for the spatial analysis of geological fractures using shapefiles.',
      keywords='geological spatial fracture analysis',
      url='https://github.com/TerminusEst/geo_frac_analysis',
      author='Sean Blake',
      author_email='blakese@tcd.ie',
      license='MIT',
      packages=['geo_frac_analysis'],
      install_requires=[
            'pyshp', 'numpy', 'matplotlib'
            ],
      include_package_data=True,
      zip_safe=False)

