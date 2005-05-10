from distutils.core import setup, Extension

compiler="linux32_nag_gcc"
if compiler=="mac_nag_gcc":
    numArrayIncludePath="/darwinports/include/python/numarray"
    extra_link_args = ['-v','/usr/local/lib/NAGWare/quickfit.o',
                                       '/usr/local/lib/NAGWare/libf97.dylib',
                                       '/usr/local/lib/NAGWare/libf96.a', '-lm',
                                       '-framework','vecLib',
                                       '-Xlinker','-Y','-Xlinker','10']
elif compiler=="linux32_nag_gcc":
    numArrayIncludePath="/home/fawzi/include/python/numarray"
    extra_link_args = [
        "-v","/apps/NAGWare/lib/quickfit.o",
        "-L/home/fawzi/cp2k/lib/Linux-i686-nag/sdbg",
        "/home/fawzi/cp2k/lib/Linux-i686-nag/sdbg/libcp2k_lib.a",
        "-lcp2k_base_lib", "-llapack", "-lg2c", "-Wl,-rpath,/apps/NAGWare/lib",
        "/apps/NAGWare/lib/libf96.so", "/apps/NAGWare/lib/libf96.a", "-lm"]
    
module1 = Extension('cp2k_interface_low',
                    include_dirs = ['/darwinports/include/python2.3'
                                    ,numArrayIncludePath,'..'
                                    ],
                    libraries = [],
                    library_dirs = [],
                    runtime_library_dirs = [],
                    extra_objects = [],
                    extra_compile_args = ['-g','-Wall','-pedantic',
                                          '-DCOMP_LINUX32_NAG_GCC'],
                    extra_link_args = extra_link_args,
                    sources = ["cp2k_c_bridge.c","cp2k_interface.c"])


setup (name = 'cp2k_interface_low',
       version = '1.0',
       description = 'This is the low level interface for CP2k',
       ext_modules = [module1])
