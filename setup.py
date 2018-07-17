from distutils.core import setup

setup(
    name='MiCall', 
    packages=['micall', 'micall.utils'],
    package_data={'micall': ['projects.json', 'project_scoring.json']}
)
