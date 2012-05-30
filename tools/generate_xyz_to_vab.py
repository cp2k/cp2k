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
#  with less agressive loop unrolling can be specified
#
# Usage:
#   generate_xyz_to_vab <template> <destination> <optional>
#
# The range of values can be specified
#
######################################


# Range of values to generate the routine variants
la_max_local_init = 0
la_max_local_end = 4
lb_max_local_init = 0
lb_max_local_end = 4

# ~~~~~~~~~~~~~~~~
# Import the template subsystem
import pyratemp

import sys

def create_routine_template(template_filename):
        """ Load the routine template string from a file and return the object """
	template_str = "\n".join([e.rstrip() for e in open(template_filename).readlines()])
	return  pyratemp.Template(template_str)

def create_caller_template():
        """ Load the caller routine template and return the object """
	call_template_str ="""
	    SUBROUTINE call_to_xyz_to_vab

	<$--(if cray_tracing)-->
	    include "pat_apif.h"
	    INTEGER :: istat
	    CALL PAT_sampling_state(PAT_STATE_ON,istat)
	    CALL PAT_region_begin(1,"XYZ_TO_VAB_X_X",istat)
	<$--(end)-->   

	    IF (la_max_local > @<la_max_local>@  .OR. lb_max_local > @<lb_max_local>@) THEN
		CALL xyz_to_vab
	    <$--(for la_value,lb_value in pairs)-->
	      ELSEIF ((la_max_local == @<la_value>@) .AND. (lb_max_local == @<lb_value>@)) THEN
		CALL xyz_to_vab_@<la_value>@_@<lb_value>@
	    <$--(end)-->
	      ENDIF

	<$--(if cray_tracing)-->
	    CALL PAT_region_end(1, istat)
	    CALL PAT_sampling_state(PAT_STATE_OFF,istat)
	<$--(end)-->

	    END SUBROUTINE call_to_xyz_to_vab
	"""
	return pyratemp.Template(call_template_str)


def create_caller_template():
        """ Load the caller routine template and return the object """
	call_template_str ="""
	    SUBROUTINE call_to_xyz_to_vab

	<$--(if cray_tracing)-->
	    include "pat_apif.h"
	    INTEGER :: istat
	    CALL PAT_sampling_state(PAT_STATE_ON,istat)
	    CALL PAT_region_begin(1,"XYZ_TO_VAB_X_X",istat)
	<$--(end)-->   

            SELECT CASE(la_max_local)
             <$--(for la_value in range(0,la_max_local+1))--> 
                 CASE (@<la_value>@)
                      SELECT CASE(lb_max_local)
                      <$--(for lb_value in range(0,lb_max_local+1))--> 
                               CASE(@<lb_value>@) 
                                    CALL xyz_to_vab_@<la_value>@_@<lb_value>@
                      <$--(end)-->  
                               CASE DEFAULT
                                    CALL xyz_to_vab 
                       END SELECT
             <$--(end)-->  
                 CASE DEFAULT
                    CALL xyz_to_vab 
             END SELECT
	   
	<$--(if cray_tracing)-->
	    CALL PAT_region_end(1, istat)
	    CALL PAT_sampling_state(PAT_STATE_OFF,istat)
	<$--(end)-->

	    END SUBROUTINE call_to_xyz_to_vab
	"""
	return pyratemp.Template(call_template_str)

def create_ordered_pairs():
        """ Create an ordered list with the most common combinations to the
            XYZ_TO_VAB routine at the beggining.
	    Most common values are:
		(1, 0)--4275300
		(1, 1)--4217600
		(0, 0)--3878400
		(0, 1)--3750300
        """
	# Order the most common cases
	l = [ (1,0) , (1,1), (0,0), (0,1) ]
        # The rest are appended as is
	for i in range(0,la_max_local_end+1):
	    for j in range(0,lb_max_local_end+1):
		if not (i,j) in l:
		    l.append((i,j,))
	return l

instance_header = """
! *****************************************************************************
!> \\brief provides specific implementations for common cases of xyz_to_vab routine.
!> 
!> \\note
!>     ____              _ _     __  __           _ _  __         _____ _     _       _____ _ _      _
!>    |  _ \  ___  _ __ ( ) |_  |  \/  | ___   __| (_)/ _|_   _  |_   _| |__ (_)___  |  ___(_) | ___| |
!>    | | | |/ _ \| '_ \|/| __| | |\/| |/ _ \ / _` | | |_| | | |   | | | '_ \| / __| | |_  | | |/ _ \ |
!>    | |_| | (_) | | | | | |_  | |  | | (_) | (_| | |  _| |_| |   | | | | | | \__ \ |  _| | | |  __/_|
!>    |____/ \___/|_| |_|  \__| |_|  |_|\___/ \__,_|_|_|  \__, |   |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>                                                        |___/
!>      ____ _                  ___                              _ _       _       _
!>     / ___| | ___  ___  ___  |_ _|_ __ ___  _ __ ___   ___  __| (_) __ _| |_ ___| |_   _
!>    | |   | |/ _ \/ __|/ _ \  | || '_ ` _ \| '_ ` _ \ / _ \/ _` | |/ _` | __/ _ \ | | | |
!>    | |___| | (_) \__ \  __/  | || | | | | | | | | | |  __/ (_| | | (_| | ||  __/ | |_| |
!>     \____|_|\___/|___/\___| |___|_| |_| |_|_| |_| |_|\___|\__,_|_|\__,_|\__\___|_|\__, |
!>                                                                                   |___/
!>     _____ _     _       _____ _ _      _
!>    |_   _| |__ (_)___  |  ___(_) | ___| |
!>      | | | '_ \| / __| | |_  | | |/ _ \ |
!>      | | | | | | \__ \ |  _| | | |  __/_|
!>      |_| |_| |_|_|___/ |_|   |_|_|\___(_)
!>
!>      This is a template
!>
!>      **** DO NOT MODIFY THIS .f90 FILE ****
!>      modify the .template instead
!>      The script to regenerate this file can be found at cp2k/tools/generate_xyz_to_vab.py
!> \par History
!>      05.2012 created
!>
!> \\author Ruyman Reyes
! *****************************************************************************
"""

if __name__ == "__main__":
	if len(sys.argv) < 3:
	    print " Incorrect number of parameters "
	    print " Usage: " + str(sys.argv[0]) + " <template_file> " + " <destination_file>"
            print " --> Using an unrolled version of the template "
	    print " Usage: " + str(sys.argv[0]) + " <template_file> " + " <destination_file> <unrolled_template_file>"
	    sys.exit(1)

	template_name = sys.argv[1]
	destination_name = sys.argv[2]
        template_nounroll = template_name
        # User can specify an unrolled version of the template
        if len(sys.argv) == 4:
            template_no_unroll = sys.argv[3]
        # Create a list of pairs ordered by probability
	l = create_ordered_pairs()
        # Initialize the templates
	call_template = create_caller_template()
	routine_template = create_routine_template(template_name)
	routine_template_nounroll = create_routine_template(template_no_unroll)

	instance  = instance_header
	# Generate routine versions 
	for pair in l:  
                if pair[0] >= 1 and pair[1] >=1 and (pair[0]+pair[1]>4):
		    instance += routine_template_nounroll(la_max_local=pair[0],lb_max_local=pair[1])
		else:
		    instance += routine_template(la_max_local=pair[0],lb_max_local=pair[1])

	# Add the caller
	instance += call_template(cray_tracing = False, 
				 la_max_local=la_max_local_end,
				 lb_max_local=lb_max_local_end,
				 pairs = l)

	# Write everyting to the external file
	include_file = open(destination_name,'w+')
	print >>include_file, instance 
        include_file.close()
