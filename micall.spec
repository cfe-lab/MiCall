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
          console=True )
