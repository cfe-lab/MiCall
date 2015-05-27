# -*- mode: python -*-
a = Analysis(['micall.py'],
             pathex=[''],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

a.datas += [
    ('micall/projects.json', 'micall/projects.json', 'DATA'),
    ('bowtie2', 'bin/bowtie2', 'DATA'),
    ('bowtie2-build', 'bin/bowtie2-build', 'DATA')
]

a.binaries += [
    ('bowtie2-align-l', 'bin/bowtie2-align-l', 'BINARY'),
    ('bowtie2-align-s', 'bin/bowtie2-align-s', 'BINARY'),
    ('bowtie2-build-l', 'bin/bowtie2-build-l', 'BINARY'),
    ('bowtie2-build-s', 'bin/bowtie2-build-s', 'BINARY'),
    ('samtools', 'bin/samtools', 'BINARY')
]

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
          console=False )

app = BUNDLE(exe,
             name='micall.app',
             icon='gui/micall-icon.icns')
