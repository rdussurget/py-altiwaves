from distutils.core import setup

setup(
    name='py-altiwaves',
    version='0.1.0',
    author='R. Dussurget',
    author_email='renaud.dussurget@gmail.com',
    packages=['kernel', 'altimetry','altimetry.config','altimetry.data','altimetry.tools','altimetry.externals','external'],
    scripts=['bin/test_detection.py','bin/test_spectral_analysis.py'],
    url='https://code.google.com/p/py-altiwaves/',
    license='LICENSE.TXT',
    description='PyALTIWAVES: Python-based ALong-Track Inventory of WAVelet based EddieS',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy",
        "os",
        "sys",
        "netCDF4"
    ],
)