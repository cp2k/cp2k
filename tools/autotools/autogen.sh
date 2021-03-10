#!/bin/sh

rm -rf autom4te.cache
rm -f aclocal.m4 ltmain.sh

echo "Running aclocal..."
aclocal "$ACLOCAL_FLAGS" -I m4 || exit 1
echo "Running autoheader..."
autoheader || exit 1
echo "Running autoconf..."
autoconf || exit 1
echo "Running libtoolize..."
(libtoolize --copy --automake || glibtoolize --automake) || exit 1
echo "Running automake..."
automake --add-missing --copy --gnu || exit 1
