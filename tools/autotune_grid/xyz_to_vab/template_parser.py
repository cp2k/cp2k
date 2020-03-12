######################################
#
#         GENERATE_XYZ_TO_VAB
#
#
# Author: Ruyman Reyes (rreyes@epcc.ed.ac.uk)
#
# Generate unrolled version of the xyz_to_vab routine
#  for the most common la_max_local and lb_max_local combinations.
# For those combinations whose unrolled code would be larger but
#  still some optimisation is needed, an additional template
#  with less aggressive loop unrolling can be specified
#
# Usage:
#   generate_xyz_to_vab <template> <destination> <optional>
#
# The range of values can be specified
#
######################################


# Range of values to generate the routine variants
from __future__ import absolute_import
from __future__ import print_function

la_max_local_init = 0
la_max_local_end = 4
lb_max_local_init = 0
lb_max_local_end = 4

# ~~~~~~~~~~~~~~~~
# Import the template subsystem
import pyratemp

import sys


def load_template(template_filename):
    """ Load the routine template string from a file and return the object """
    # template_str = "\n".join([e.rstrip() for e in open(template_filename).readlines()])
    return pyratemp.Template(filename=template_filename)


def create_ordered_pairs():
    """ Create an ordered list with the most common combinations to the
        XYZ_TO_VAB routine at the beginning.
        Most common values are:
        (1, 0)--4275300
        (1, 1)--4217600
        (0, 0)--3878400
        (0, 1)--3750300
    """
    # Order the most common cases
    l = [(1, 0), (1, 1), (0, 0), (0, 1)]
    # The rest are appended as is
    for i in range(0, la_max_local_end + 1):
        for j in range(0, lb_max_local_end + 1):
            if not (i, j) in l:
                l.append((i, j,))
    return l


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(" Incorrect number of parameters ")
        print(
            " Usage: " + str(sys.argv[0]) + " <template_file> " + " <destination_file>"
        )
        sys.exit(1)

    template_file_name = sys.argv[1]
    output_file_name = sys.argv[2]
    # Create a list of pairs ordered by probability
    l = create_ordered_pairs()
    # Initialize the templates
    module_template = load_template(template_file_name)
    source_str = module_template(
        cray_tracing=False,
        la_max_local=la_max_local_end,
        lb_max_local=lb_max_local_end,
        pairs=l,
    )
    source_file = open(output_file_name, "w+")
    print(source_str, file=source_file)

    source_file.close()
