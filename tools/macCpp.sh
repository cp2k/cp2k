#!/bin/sh
# workaround to remove the stupid "feature" of cpp 3.3 that *always* adds a
# #pragma GCC set_debug_pwd to each file
cpp $* | sed -n '/^#pragma/!p'
