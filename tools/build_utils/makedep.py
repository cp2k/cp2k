#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re, sys
from os import path
from os.path import dirname, basename, normpath

# pre-compiled regular expressions
re_module = re.compile(r"(?:^|\n)\s*module\s+(\w+)\s.*\n\s*end\s*module",re.DOTALL)
re_program = re.compile(r"\n\s*end\s*program")
re_main   = re.compile(r"\sint\s+main\s*\(")
re_use    = re.compile(r"\n\s*use\s+(\w+)")
re_incl1  = re.compile(r"\n\s*include\s+['\"](.+)['\"]")
re_incl2  = re.compile(r"\n#include\s+['\"](.+)['\"]")


#=============================================================================
def main():
    messages = []
    #process command line arguments
    out_fn  = sys.argv[1]
    mod_format = sys.argv[2]
    mode = sys.argv[3]
    archive_ext = sys.argv[4]
    src_files = sys.argv[5:]
    if(mod_format not in ('lower', 'upper', 'no')):
        error('Module filename format must be eighter of "lower", "upper", or "no".')
    if(mode not in ('normal', 'hackdep', 'mod_compiler')):
        error('Mode must be eighter of "normal", "hackdep", or "mod_compiler".')
    for fn in src_files:
        if(not fn.startswith("/")):
            error("Path of source-file not absolut: "+fn)
    src_basenames = [basename(fn).rsplit(".", 1)[0] for fn in src_files]
    for bfn in src_basenames:
        if(src_basenames.count(bfn) > 1):
            error("Multiple source files with the same basename: "+bfn)

    # parse files
    parsed_files = dict()
    for fn in src_files:
        parse_file(parsed_files, fn) #parses also included files
    messages.append("Parsed %d files"%len(parsed_files))

    # create table mapping fortan module-names to file-name
    mod2fn = dict()
    for fn in src_files:
        for m in parsed_files[fn]['module']:
            if(mod2fn.has_key(m)):
                error('Multiple declarations of module "%s"'%m)
            mod2fn[m] = fn
    messages.append("Created mod2fn table, found %d modules."%len(mod2fn))

    # check "one module per file"-convention
    for m, fn in mod2fn.items():
        if(basename(fn) != m+".F"):
            error("Names of module and file do not match: "+fn)

    # read package manifests
    packages = dict()
    for fn in src_files:
        p = normpath(dirname(fn))
        read_pkg_manifest(packages, p)
    messages.append("Read %d package manifests"%len(packages))

    # check dependencies against package manifests
    n_deps = 0
    for fn in src_files:
        p = normpath(dirname(fn))
        if(not parsed_files[fn]['program']):
            packages[p]['objects'].append(src2obj(basename(fn)))
        deps = collect_include_deps(parsed_files, fn)
        deps += [ mod2fn[m] for m in collect_use_deps(parsed_files, fn) if mod2fn.has_key(m) ]
        n_deps += len(deps)
        for d in deps:
            dp = normpath(dirname(d))
            if(dp not in packages[p]['allowed_deps']):
                error("Dependency forbidden according to package manifest: %s -> %s"%(fn, d))
            if(dp != p and packages[dp].has_key("public")):
                if(basename(d) not in packages[dp]["public"]):
                    error("File not public according to package manifest: %s -> %s"%(fn, d))
    messages.append("Checked %d dependencies"%n_deps)

    # check for circular dependencies
    for fn in parsed_files.keys():
        find_cycles(parsed_files, mod2fn, fn)

    # write messages as comments
    makefile = "".join(["#makedep: %s\n"%m for m in messages])
    makefile += "\n"

    # write rules for archives
    for p in packages.keys():
        if(len(packages[p]['objects']) > 0):
            makefile += "# Package %s\n"%p
            makefile += "$(LIBDIR)/%s : "%(packages[p]['archive']+archive_ext)
            makefile += " ".join(packages[p]['objects']) + "\n\n"

    # write rules for executables
    archive_postfix = archive_ext.rsplit(".",1)[0]
    for fn in src_files:
        if(not parsed_files[fn]['program']):
            continue
        bfn = basename(fn).rsplit(".", 1)[0]
        makefile += "# Program %s\n"%fn
        makefile += "$(EXEDIR)/%s.$(ONEVERSION) : %s.o "%(bfn, bfn)
        p = normpath(dirname(fn))
        deps = collect_pkg_deps(packages, p)
        makefile += " ".join(["$(LIBDIR)/"+a+archive_ext for a in deps]) + "\n"
        makefile += "\t" + "$(LD) $(LDFLAGS) -L$(LIBDIR) -o $@ %s.o "%bfn
        makefile += "$(EXTERNAL_OBJECTS) "
        assert(all([a.startswith("lib") for a in deps]))
        makefile += " ".join(["-l"+a[3:]+archive_postfix for a in deps])
        makefile += " $(LIBS)\n\n"

    # write rules for objects
    for fn in src_files:
        deps = " ".join(collect_include_deps(parsed_files, fn))
        mods = collect_use_deps(parsed_files, fn)
        mods.sort(cmp_mods) # sort mods to speedup compilation
        for m in mods:
            if(mod2fn.has_key(m)):
                deps += " " + mod2modfile(m, mod_format)
        if(mode == "hackdep"):
            deps = ""

        bfn = basename(fn)
        makefile += "# Object %s\n"%bfn
        provides = [mod2modfile(m, mod_format) for m in parsed_files[fn]['module']]
        for mfn in provides:
            makefile += "%s : %s "%(mfn, bfn) + deps + "\n"
        makefile += "%s : %s "%(src2obj(bfn), bfn) + deps
        if(mode == "mod_compiler"):
            makefile += " " + " ".join(provides)
        makefile += "\n\n"

    f = open(out_fn, "w")
    f.write(makefile)
    f.close()


#=============================================================================
def cmp_mods(a,b):
    # list "type" modules first, they are probably on the critical path
    if("type" in a):
        return -1
    if("type" in b):
        return +1
    return 0


#=============================================================================
def parse_file(parsed_files, fn):
    if(fn in parsed_files): return

    content = open(fn).read()

    # re.IGNORECASE is horribly expensive. Converting to lower-case upfront
    content_lower = content.lower()

    # all files are parsed for cpp includes
    incls = re_incl2.findall(content) #cpp includes, those are case-sensitiv

    mods=[]; uses=[]; prog=False;
    if(fn[-2:]==".F" or fn[-4:]==".f90"):
        mods += re_module.findall(content_lower)
        prog  = re_program.search(content_lower) != None
        uses += re_use.findall(content_lower)
        incl1_iter = re_incl1.finditer(content_lower) # fortran includes
        incls += [ content[m.start(1):m.end(1)] for m in incl1_iter ]

    if(fn[-2:] == ".c" or fn[-3:]==".cu"):
        prog = re_main.search(content) != None # C is case-sensitiv

    # exclude included files from outside the source tree
    def incl_fn(i):
        return normpath(path.join(dirname(fn), i))

    existing_incl = [i for i in incls if path.exists(incl_fn(i))]

    # store everything in parsed_files cache
    parsed_files[fn] = {'module':mods, 'program': prog, 'use':uses, 'include':existing_incl}

    # parse included files
    for i in existing_incl:
        parse_file(parsed_files, incl_fn(i))


#=============================================================================
def read_pkg_manifest(packages, p):
    if(packages.has_key(p)): return

    fn = p+"/PACKAGE"
    if(not path.exists(fn)):
        error("Could not open PACKAGE manifest: "+fn)
    content = open(fn).read()

    packages[p] = eval(content)
    packages[p]['objects'] = []
    if(not packages[p].has_key("archive")):
        packages[p]['archive'] = "libcp2k"+basename(p)
    packages[p]['allowed_deps'] = [normpath(p)]
    packages[p]['allowed_deps'] += [normpath(path.join(p,r)) for r in packages[p]['requires']]

    for r in packages[p]['requires']:
        read_pkg_manifest(packages, normpath(path.join(p,r)))


#=============================================================================
def mod2modfile(m, mod_format):
    if(mod_format == 'no'):
        return("")
    if(mod_format == 'lower'):
        return(m.lower() + ".mod")
    if(mod_format == 'upper'):
        return(m.upper() + ".mod")
    assert(False) # modeps unknown


#=============================================================================
def src2obj(src_fn):
    return( basename(src_fn).rsplit(".",1)[0] + ".o" )


#=============================================================================
def collect_include_deps(parsed_files, fn):
    pf = parsed_files[fn]
    incs = []

    for i in pf['include']:
        fn_inc = normpath(path.join(dirname(fn), i))
        if(parsed_files.has_key(fn_inc)):
            incs.append(fn_inc)
            incs += collect_include_deps(parsed_files, fn_inc)

    return(list(set(incs)))


#=============================================================================
def collect_use_deps(parsed_files, fn):
    pf = parsed_files[fn]
    uses = pf['use']

    for i in pf['include']:
        fn_inc = normpath(path.join(dirname(fn), i))
        if(parsed_files.has_key(fn_inc)):
            uses += collect_use_deps(parsed_files, fn_inc)

    return(list(set(uses)))


#=============================================================================
def find_cycles(parsed_files, mod2fn, fn, S=None):
    pf = parsed_files[fn]
    if(pf.has_key('visited')):
        return

    if(S==None):
        S = []

    for m in pf['module']:
        if(m in S):
            i = S.index(m)
            error("Circular dependency: "+ " -> ".join(S[i:] + [m]))
        S.append(m)

    for m in collect_use_deps(parsed_files, fn):
        if(mod2fn.has_key(m)):
            find_cycles(parsed_files, mod2fn, mod2fn[m], S)

    for m in pf['module']:
        S.pop()

    pf['visited'] = True


#=============================================================================
def collect_pkg_deps(packages, p, archives=None, S=None):
    if(archives == None): archives = []
    if(S == None): S = []

    a = packages[p]['archive']
    if(a in archives):
        return(archives)

    if(a in S):
        i = S.index(a)
        error("Circular package dependency: "+ " -> ".join(S[i:] + [a]))
    S.append(a)

    for r in packages[p]['requires']:
        d = normpath(path.join(p,r))
        collect_pkg_deps(packages, d, archives, S)

    S.pop()

    if(len(packages[p]['objects']) > 0):
        archives.insert(0, packages[p]['archive'])

    return(archives)


#=============================================================================
def error(msg):
    sys.stderr.write("makedep error: %s\n"%msg)
    sys.exit(1)

#=============================================================================
# Python 2.4 compatibility
def all(iterable):
    for element in iterable:
        if(not element):
            return(False)
    return(True)

#=============================================================================
main()

#EOF
