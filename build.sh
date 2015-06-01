#!/bin/sh
# use PyInstaller to bundle MiCall into an OS-X .app executable
pyinstaller -F -w --icon=gui/grabber.icns --clean micall.spec --noconfirm
