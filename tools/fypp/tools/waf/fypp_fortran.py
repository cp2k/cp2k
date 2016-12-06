#!/usr/bin/env python
# encoding: utf-8
# Bálint Aradi, 2016

'''Uses Fypp as Fortran preprocessor (.fpp -> .f90). Use this one (instead of 
fypp_preprocessor) if you want to preprocess Fortran sources with Fypp.

You typically trigger the preprocessing via the 'fypp' feature::

	def options(opt):
		opt.load('compiler_c')
		opt.load('compiler_fc')
		opt.load('fypp_fortran')
	
	def configure(conf):
		conf.load('compiler_c')
		conf.load('compiler_fc')
		conf.load('fypp_fortran')
	
	
	def build(bld):
		sources = bld.path.ant_glob('*.fpp')
		
		bld(
			features='fypp fc fcprogram',
			source=sources,
			target='myprog'
		)

Please check the documentation in the fypp_preprocessor module for the
description of the uselib variables which may be passed to the task generator.
'''

from waflib import TaskGen
import fypp_preprocessor


################################################################################
# Configure
################################################################################

def configure(conf):
	fypp_preprocessor.configure(conf)
	

################################################################################
# Build
################################################################################

class fypp_fortran(fypp_preprocessor.fypp_preprocessor):
		
	ext_in = [ '.fpp' ]
	ext_out = [ '.f90' ]


@TaskGen.extension('.fpp')
def fypp_preprocess_fpp(self, node):
	'Preprocess the .fpp files with Fypp.'

	f90node = node.change_ext('.f90')
	self.create_task('fypp_fortran', node, [ f90node ])
	if 'fc' in self.features:
		self.source.append(f90node)
