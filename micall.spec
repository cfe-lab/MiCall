# -*- mode: python -*-
a = Analysis(['micall.py'],
             pathex=[''],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

a.datas += [('micall/projects.json', 'micall/projects.json', 'DATA'),
    ('bowtie2', 'bowtie2/bowtie2', 'DATA'),
    ('bowtie2-build', 'bowtie2/bowtie2-build', 'DATA')
]

a.binaries += [('bowtie2-align-l', 'bowtie2/bowtie2-align-l', 'BINARY'),
               ('bowtie2-align-s', 'bowtie2/bowtie2-align-s', 'BINARY'),
               ('bowtie2-build-l', 'bowtie2/bowtie2-build-l', 'BINARY'),
               ('bowtie2-build-s', 'bowtie2/bowtie2-build-s', 'BINARY')]

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='micall',
          debug=False,
          strip=None,
          upx=True,
          console=True )
#app = BUNDLE(exe,
#             name='micall.app',
#             icon='micall/gui/grabber.icns')
