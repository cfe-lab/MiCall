#! /bin/sh

CURDIR="${0%/*}"

set -xe

cp -- "$CURDIR"/fakeiva.py /bin/iva
chmod +x -- /bin/iva

mkdir -- /data
cp -T -r -- "$CURDIR"/build/assembly-ios /data/assembly-ios
