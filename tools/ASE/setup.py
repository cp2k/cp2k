from distutils.core import setup, Extension

module1 = Extension('cp2k_interface',
                    include_dirs = ['/darwinports/include/python2.3'
                                    ,'/darwinports/include/python2.3/numarray','..'
                                    ],
                    libraries = ['cp2k_lib', 'cp2k_base_lib'],
                    library_dirs = ['/Users/fawzi/cp2k/lib/Darwin-PowerMacintosh-nag/sdbg','/usr/local/lib/NAGWare'],
                    runtime_library_dirs = ['/Users/fawzi/cp2k/lib/Darwin-PowerMacintosh-nag/sdbg','/usr/local/lib/NAGWare'],
                    extra_objects = [],
                    extra_compile_args = ['-g','-Wall','-pedantic'],
                    extra_link_args = ['-v','/usr/local/lib/NAGWare/quickfit.o',
                                       '/usr/local/lib/NAGWare/libf97.dylib',
                                       '/usr/local/lib/NAGWare/libf96.a', '-lm',
                                       '-framework','vecLib',
                                       '-Xlinker','-Y','-Xlinker','10'],
                    sources = ['cp2k_interface.c','cp2k_c_bridge.c'])


setup (name = 'cp2k_interface',
       version = '1.0',
       description = 'This is the interface for CP2k',
       ext_modules = [module1])
