# -*- mode: python -*-
a = Analysis(['micall.py'],
             pathex=[''],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

# http://stackoverflow.com/a/20695056/4794
a.datas = list({tuple(map(str.upper, t)) for t in a.datas})

a.datas += [
    ('micall/projects.json', 'micall/projects.json', 'DATA'),
    ('micall/g2p/g2p.matrix', 'micall/g2p/g2p.matrix', 'DATA'),
    ('micall/g2p/g2p_fpr.txt', 'micall/g2p/g2p_fpr.txt', 'DATA'),
]

a.binaries += [
    ('bin/bowtie2-align-s.exe', 'bin/bowtie2-align-s.exe', 'BINARY'),
    ('bin/bowtie2-build-s.exe', 'bin/bowtie2-build-s.exe', 'BINARY')]

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='micall.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False,
          icon='gui/micall.ico')
