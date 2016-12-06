#!/usr/bin/env python3
import sys
import re
import os

VERSION_PATTERN = r'\d+\.\d+(?:\.\d+)?'
FILES_PATTERNS = [ ('bin/fypp', 
                    r'^VERSION\s*=\s*([\'"]){}\1'.format(VERSION_PATTERN), 
                    "VERSION = '{version}'"),
                   ('docs/fypp.rst',
                    r'Fypp Version[ ]*{}.'.format(VERSION_PATTERN),
                    'Fypp Version {shortversion}.'),
                   ('setup.py',
                    r'version\s*=\s*([\'"]){}\1'.format(VERSION_PATTERN),
                    "version='{version}'"),
                   ('docs/conf.py',
                    r'version\s*=\s*([\'"]){}\1'.format(VERSION_PATTERN),
                    "version = '{shortversion}'"),
                   ('docs/conf.py',
                    r'release\s*=\s*([\'"]){}\1'.format(VERSION_PATTERN),
                    "release = '{version}'"),
]

if len(sys.argv) < 2:
    print("Missing version string")
    sys.exit(1)


version = sys.argv[1]
shortversion = '.'.join(version.split('.')[0:2])

match = re.match(VERSION_PATTERN, version)
if match is None:
    print("Invalid version string")
    sys.exit(1)

rootdir = os.path.join(os.path.dirname(sys.argv[0]), '..')
for fname, regexp, repl in FILES_PATTERNS:
    fname = os.path.join(rootdir, fname)
    print("Replacments in '{}': ".format(fname), end='')
    fp = open(fname, 'r')
    txt = fp.read()
    fp.close()
    replacement = repl.format(version=version, shortversion=shortversion)
    newtxt, nsub = re.subn(regexp, replacement, txt, flags=re.MULTILINE)
    print(nsub)
    fp = open(fname, 'w')
    fp.write(newtxt)
    fp.close()

    
# Replace version number in Change Log and adapt decoration below
fname = os.path.join(rootdir, 'CHANGELOG.rst')
print("Replacments in '{}': ".format(fname), end='')
fp = open(fname, 'r')
txt = fp.read()
fp.close()
decoration = '=' * len(version)
newtxt, nsub = re.subn('^Unreleased\s*\n=+', version + '\n' + decoration, txt,
                       flags=re.MULTILINE)
print(nsub)
fp = open(fname, 'w')
fp.write(newtxt)
fp.close()
