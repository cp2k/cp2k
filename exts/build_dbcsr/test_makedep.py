#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unittests for the makedep.py script
"""

import unittest
import tempfile
import shutil
import json
import io
import os
from os import path

import makedep


def _mkpkg(content, base_dir):
    """
    Create a PACKAGE file in the given directory.
    If the content is a list or a dict, serialize the structure using json.dump.
    """

    with open(path.join(base_dir, "PACKAGE"), "w") as fhandle:
        if isinstance(content, (dict, list)):
            json.dump(content, fhandle, ensure_ascii=False)
        else:
            fhandle.write(content)


class TestCheckArchives(unittest.TestCase):
    """Test case for the main method in makedep"""

    def setUp(self):
        self.base_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.base_dir)

    def test_empty_dir(self):
        """
        running on an empty dir with no files generates an empty deps file
        """
        my_dir = path.join(self.base_dir, "empty")
        out_fn = path.join(my_dir, "all.dep")

        os.mkdir(my_dir)

        makedep.main(out_fn, "empty_project", "lower", "normal", ".a", my_dir, [])

        with open(out_fn, "r") as fhandle:
            no_comment_lines = [
                l.strip() for l in fhandle if not l.startswith("#") and l.strip()
            ]

        self.assertEqual(len(no_comment_lines), 0)

    def test_single_empty(self):
        """
        running on a dir with a single file
        """

        my_dir = path.join(self.base_dir, "single_empty")
        out_fn = path.join(my_dir, "all.dep")
        single_fn = path.join(my_dir, "single.F")

        os.mkdir(my_dir)
        with open(single_fn, "w") as fhandle:
            fhandle.write("\n")

        with self.assertRaises(SystemExit):
            # should throw an exception due to missing PACKAGES
            makedep.main(
                out_fn,
                "faulty_project",
                "lower",
                "normal",
                ".a",
                my_dir,
                ["./single.F"],
            )

        _mkpkg(
            {
                "description": "Nothing",
                "archive": "test",
                "public": ["*.F"],
                "requires": [],
            },
            my_dir,
        )

        makedep.main(
            out_fn, "single_empty", "lower", "normal", ".a", my_dir, ["./single.F"]
        )

        with open(out_fn, "r") as fhandle:
            no_comment_lines = [
                l.strip() for l in fhandle if not l.startswith("#") and l.strip()
            ]

        self.assertEqual(
            no_comment_lines,
            [
                "$(LIBDIR)/test.a : single.o",
                "install: PUBLICFILES += *.F",
                "single.o : single.F",
            ],
        )

    def test_unicode(self):
        """
        running on a dir with a single file
        """

        my_dir = path.join(self.base_dir, "unicode")
        out_fn = path.join(my_dir, "all.dep")
        single_fn = path.join(my_dir, "single.F")

        os.mkdir(my_dir)
        with io.open(single_fn, "w", encoding="utf8") as fhandle:
            fhandle.write(u"! Ångström\n")

        _mkpkg(
            """{
    "description": "unicode test with just an Ångström 😉",
    "archive": "test",
    "public": ["*.F"],
    "requires": []
}""",
            my_dir,
        )

        makedep.main(
            out_fn, "unicode_project", "lower", "normal", ".a", my_dir, ["./single.F"]
        )

        with open(out_fn, "r") as fhandle:
            no_comment_lines = [
                l.strip() for l in fhandle if not l.startswith("#") and l.strip()
            ]

        self.assertEqual(
            no_comment_lines,
            [
                "$(LIBDIR)/test.a : single.o",
                "install: PUBLICFILES += *.F",
                "single.o : single.F",
            ],
        )

    def test_subpackage(self):
        """
        running on a dir with a single file
        """

        my_dir = path.join(self.base_dir, "subpackage")
        out_fn = path.join(my_dir, "all.dep")
        sub_dir = path.join(my_dir, "sub")
        single_fn = path.join(my_dir, "single.F")
        sub_fn = path.join(sub_dir, "sub.F")

        os.mkdir(my_dir)
        os.mkdir(sub_dir)

        for fname in [single_fn, sub_fn]:
            with open(fname, "w") as fhandle:
                fhandle.write("")

        _mkpkg(
            """{
    "description": "test with a simple subpackage",
    "archive": "test",
    "public": ["*.F"],
    "requires": ["./sub"]
}""",
            my_dir,
        )

        _mkpkg(
            """{
    "description": "the sub package",
    "requires": [""]
}""",
            sub_dir,
        )

        makedep.main(
            out_fn,
            "sub_pkg_without_dep",
            "lower",
            "normal",
            ".a",
            my_dir,
            ["./single.F"],
        )

        with open(out_fn, "r") as fhandle:
            no_comment_lines = [
                l.strip() for l in fhandle if not l.startswith("#") and l.strip()
            ]

        self.assertEqual(
            no_comment_lines,
            [
                "$(LIBDIR)/test.a : single.o",
                "install: PUBLICFILES += *.F",
                "single.o : single.F",
            ],
        )

        with open(single_fn, "w") as fhandle:
            fhandle.write(
                """module single
    use :: sub
end module"""
            )

        with open(sub_fn, "w") as fhandle:
            fhandle.write(
                """module sub
end module"""
            )

        makedep.main(
            out_fn,
            "sub_pkg_with_dep",
            "lower",
            "normal",
            ".a",
            my_dir,
            ["./single.F", "./sub/sub.F"],
        )

        with open(out_fn, "r") as fhandle:
            no_comment_lines = [
                l.strip() for l in fhandle if not l.startswith("#") and l.strip()
            ]

        self.assertEqual(
            sorted(no_comment_lines),
            [
                "$(LIBDIR)/libsub_pkg_with_depsub.a : sub.o",
                "$(LIBDIR)/test.a : single.o",
                "install: PUBLICFILES += *.F",
                "single.mod : single.F",
                "single.o : single.F",
                "sub.o : sub.F",
            ],
        )


# ============================================================================
if __name__ == "__main__":
    unittest.main()
