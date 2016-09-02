from distutils.core import setup

setup(
    name='MiCall', 
    packages=['micall', 'micall.core', 'micall.utils'],
    package_data={'micall': ['projects.json', 'project_scoring.json']}
)
