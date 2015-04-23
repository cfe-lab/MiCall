# -*- mode: python -*-
a = Analysis(['micall.py'],
             pathex=['/Users/art/git/MiCall'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas+[('micall/projects.json', 'micall/projects.json', 'DATA')],
          name='micall',
          debug=False,
          strip=None,
          upx=True,
          console=True )
app = BUNDLE(exe,
             name='micall.app',
             icon=None)
