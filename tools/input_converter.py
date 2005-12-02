#!env python
import sys, os, os.path, StringIO
from input_parser import *

def guaranteePath(root_sect,path):
    p=path.split("%")
    s_path=[root_sect]
    for s in p:
        if s_path[-1].subsections.has_key(s.lower()):
            s_path.append(s_path[-1].subsections[s.lower()][0])
        else:
            new_sect=Section(s)
            new_sect.auto_g=1
            s_path[-1].add_subsection(new_sect)
            s_path.append(new_sect)
    return s_path

class Conv:
    def __init__(self,old_name,path,upcase=str.upper,sort_keywords=1,special_convert=None,
                 merge=1):
        self.old_name=old_name
        self.path=path
        self.upcase=upcase
        self.sort_keywords=sort_keywords
        self.special_convert=special_convert
        self.merge=merge

class SectionConverter:
    def __init__(self):
        self.conversions={}
    def add_conversion(self,conv):
        if self.conversions.has_key(conv.old_name.lower()):
            raise KeyError("conversion for "+conv.old_name+" already present")
        self.conversions[conv.old_name.lower()]=conv
    def convert(self,oldInput,newInput):
        newInput.pre_comments=oldInput.pre_comments
        newInput.post_comments=oldInput.post_comments
        for sects in oldInput.subsections.itervalues():
            for sect in sects:
                if not self.conversions.has_key(sect.name.lower()):
                    print "WARNING dropping section",sect.name,"unknown conversion"
                else:
                    conv=self.conversions[sect.name.lower()]
                    sect_att=None
                    if conv.path:
                        path=conv.upcase(conv.path)
                        s_att=guaranteePath(newInput,path)
                        if conv.merge:
                            sect_att=self.merge_section(sect,s_att[-1],
                                                        s_att[-2],conv)
                        else:
                            sect_att=s_att[-1]
                    if conv.special_convert:
                        conv.special_convert(oldInput,sect,newInput,sect_att,conv)
        self.globalCleanup(newInput)
    
    def merge_section(self,sect,sect_att,super_sect,conv):
        if sect_att and hasattr(sect_att,"auto_g"):
            del sect_att.auto_g
        else:
            new_sect=Section(sect.name)
            super_sect.add_subsection(new_sect)
            sect_att=new_sect
        k_list=sect.keywords_list
        for kw in sect.keywords_list:
            kw.name=conv.upcase(kw.name)
            sect_att.add_keyword(kw)
        if conv.sort_keywords:
            sect_att.sort_keywords()
        for subSs in sect.subsections.values():
            for subS in subSs:
                subS.name=conv.upcase(subS.name)
                if sect_att.subsections.has_key(subS.name.lower()):
                    sub_sect=sect_att.subsections[subS.name.lower()][0]
                    if not hasattr(sub_sect,'auto_g'):
                        if subS.keywords_list:
                            print "WARNING duplicate non firstlevel section",subS.name
                            print "        having keywords, mixing new/old input?"
                        else:
                            print "NOTE fusing keywordless section",subS.name
                            sub_sect.auto_g=1
                else:
                    sub_sect=None
                self.merge_section(subS,sub_sect,sect_att,conv)
        sect_att.args+=sect.args
        sect_att.pre_comments+=sect.pre_comments
        sect_att.post_comments+=sect.post_comments
        return sect_att

    def globalCleanup(self,newInput):
        return None

mainConv=SectionConverter()

def dft_conv(oldInput,oldSect,newInput,new_sect,conv):
    toXc=["density_cutoff","gradient_cutoff","tau_cutoff",
          "density_smooth_cutoff_range","functional_routine"]
    for k in new_sect.keywords_list:
        kname=k.name.lower()
        if kname[:3]=="sic":
            s=guaranteePath(new_sect,conv.upcase("sic"))
            if kname=="sic":
                s[-1].add_keyword(Keyword("sic_method",[k.values[0]]))
            else:
                s[-1].add_keyword(k)
            del new_sect.keywords[kname]
        elif kname=='functional':
            p=guaranteePath(newInput,conv.upcase("force_eval%dft%xc%xc_functional"))
            del p[-1].auto_g
            p[-1].args=new_sect.keywords['functional'][0].values
            if p[-1].args[0].lower()=="my_tpss": p[-1].args[0]="tpss"
            del new_sect.keywords['functional']
        elif kname=="xc_deriv" or kname=="xc_smooth_rho":
            s=guaranteePath(new_sect,conv.upcase("xc%xc_grid"))
            s[-1].add_keyword(k)
            del new_sect.keywords[kname]
        elif (kname in toXc):
            s=guaranteePath(new_sect,conv.upcase("xc"))
            s[-1].add_keyword(k)
            del new_sect.keywords[kname]
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("DFT","FORCE_EVAL%DFT",special_convert=dft_conv))

mainConv.add_conversion(Conv("FORCE_EVAL","FORCE_EVAL"))

mainConv.add_conversion(Conv("METADYNAMICS","MOTION%METADYN"))
mainConv.add_conversion(Conv("COLVAR","MOTION%METADYN%COLVAR"))

def cp2k_conv(oldInput,oldSect,newInput,new_sect,conv):
    if new_sect.keywords.has_key('iolevel'):
        del new_sect.keywords['iolevel']
        new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("CP2K","GLOBAL",special_convert=cp2k_conv))

def md_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    convert={'__default__':{'__section_path__':"motion%md",'__cut__':5,
                            'ensem':'ensemble','steps':'steps','tempe':'temperature',
                            'resta':'restart','temp_':'temp_tol','times':'timestep',
                            'gamma':'gamma','const':None},
             'nose_para':{'__section_path__':'motion%md%nose','__cut__':5,
                          'lengt':'length','yoshi':'yoshida','timec':'timecon','mts:':'mts'},
             'uniaxial_':{'__section_path__':'motion%md%uniaxial','__cut__':5,
                          'cmass':'cmass','v_sho':'v_shock','press':'pressure'},
             'hmc':{'__section_path__':"motion%md%hmc",'__cut__':80,
                    'ext_temp':'ext_temp', 'ld_steps':'ld_steps', 'qs_accept':'qs_accept',
                    'mc_method':'mc_method','md_method':'md_method','qs_start':'qs_start',
                    'recover_momenta':'recover_momenta', 'restore_history':'restore_history',
                    'rnd_velocities':'rnd_velocities', 'semi_hybrid':'semi_hybrid',
                    'temp_alpha':'temp_alpha', 'temp_beta':'temp_beta',
                    'tot_energy':'tot_energy','vel_scaling':'vel_scaling'},
             'barostat_':{'__section_path__':'motion%md%barostat','__cut__':5,
                          'press':'pressure','timec':'timecon'},
             'shake:':{'__section_path__':"force_eval%subsys%constraint",'__cut__':5,
                            'shake':'shake'}
             }
    first_check=convert['__default__']
    
    def add_keys(c,k):
        if c.has_key(kname[:c['__cut__']]):
            s=guaranteePath(newInput,conv.upcase(c['__section_path__']))
            if hasattr(s[-1],'auto_g'):
                del s[-1].auto_g
            if c[kname[:c['__cut__']]]:
                s[-1].add_keyword(Keyword(conv.upcase(c[kname[:c['__cut__']]]),
                                          values=[conv.upcase(k.values[0])]))
            if len(k.values)<2:
                first_check=convert['__default__']
            return 1
        else:
            return 0

    for k in oldSect.keywords_list:
        kname=k.name.lower()

        if not add_keys(first_check,k):
            if convert.has_key(kname[:9]):
                add_keys(convert[kname[:9]],k)
                first_check=convert[kname[:9]]
            elif not add_keys(convert['__default__'],k):
                print "WARNING ignoring unknown keyword",k.name,"in section MD"
mainConv.add_conversion(Conv("MD","MOTION%MD",special_convert=md_conv,merge=0))

def cons_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    new_sect.pre_comments=oldSect.pre_comments
    new_sect.post_comments=oldSect.post_comments
    consTypes={'g3x3':{'atoms':[3,4,5],'distances':[6,7,8],'molecule':[2]},
               'g4x6':{'atoms':[3,4,5,6],'distances':[7,8,9,10,11,12],'molecule':[2]},
               'dist':{'atoms':[3,4],'distance':[5],'molecule':[2]}}
    for line in oldSect.raw_lines:
        if line.isspace():continue
        ll=line.split()
        sname=ll[0].lower()
        if sname.lower()=='dist':
            s=guaranteePath(new_sect,conv.upcase("internals"))
            sAtt=Section(conv.upcase('distance'))
            s[-1].add_subsection(sAtt)
        else:
            sAtt=Section(conv.upcase(sname))
            new_sect.add_subsection(sAtt)
        convAtt=consTypes[sname]
        for (k,v) in convAtt.iteritems():
            kw=Keyword(conv.upcase(k),map(lambda x:ll[x],v))
            sAtt.add_keyword(kw)
mainConv.add_conversion(Conv("CONSTRAINTS","FORCE_EVAL%SUBSYS%CONSTRAINT",
                             special_convert=cons_conv,merge=0))
    
def qs_conv(oldInput,oldSect,newInput,new_sect,conv):
    s=guaranteePath(newInput,conv.upcase("FORCE_EVAL%DFT%MGRID"))
    toMg=['cutoff','ngrids','ngrid','progression_factor',
          'commensurate','realspace','rel_cutoff']
    for k in new_sect.keywords_list:
        if k.name.lower() in toMg:
            del new_sect.keywords[k.name.lower()]
            if k.name.lower()=='ngrid': k.name=conv.upcase('ngrids')
            s[-1].add_keyword(k)
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("QS","FORCE_EVAL%DFT%QS",special_convert=qs_conv))

def scf_conv(oldInput,oldSect,newInput,new_sect,conv):
    for k in new_sect.keywords_list:
        if k.name.lower()=="cholesky_off":
            del new_sect.keywords["cholesky_off"]
            if new_sect.keywords.has_key("cholesky"):
                new_sect.keywords["cholesky"][0].values=[conv.upcase('off')]
            else:
                new_sect.add_keyword(Keyword(conv.upcase("cholesky"),
                                             values=[conv.upcase('off')]))
        if k.name.lower()=="guess":
            k.name=conv.upcase("scf_guess")
            new_sect.keywords[k.name.lower()]=[k]
            del new_sect.keywords["guess"]
        elif k.name.lower()=="ot":
            s=guaranteePath(new_sect,conv.upcase("ot"))
            s[-1].args=["ON"]
            del new_sect.keywords["ot"]
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("SCF","FORCE_EVAL%DFT%SCF",
                             special_convert=scf_conv))

def print_conv(oldInput,oldSect,newInput,new_sect,conv):
    if oldSect.args:
        s=guaranteePath(newInput,conv.upcase("global"))
        if s[-1].keywords.has_key("iolevel"):
            del s[-1].keywords["iolevel"]
            s[-1].rebuild_keyword_list()
        if s[-1].keywords.has_key("print_level"):
            s[-1].keywords["print_level"][0].values=[conv.upcase(oldSect.args[0])]
        else:
            s[-1].add_keyword(Keyword(conv.upcase("print_level"),
                                      [conv.upcase(oldSect.args[0])]))
    for k in oldSect.keywords_list:
        print "WARNING dropping printkey",k.name
mainConv.add_conversion(Conv("PRINT",None,special_convert=print_conv))

mainConv.add_conversion(Conv("KIND","FORCE_EVAL%SUBSYS%KIND"))

mainConv.add_conversion(Conv("CELL","FORCE_EVAL%SUBSYS%CELL"))

mainConv.add_conversion(Conv("FARMING","FARMING"))

def jobs_conv(oldInput,oldSect,newInput,new_sect,conv):
    s=guaranteePath(newInput,conv.upcase("farming%job"))
    for line in oldSect.raw_lines:
        if line:
            l=line.split()
            if len(l)==2:
                if s[-1].auto_g:
                    del s[-1].auto_g
                else:
                    s[-1]=Section("JOB")
                    s[-2].add_subsection(s[-1])
                s[-1].add_keyword(Keyword(conv.upcase("directory"),
                                          [l[0]]))
                s[-1].add_keyword(Keyword(conv.upcase("input_file_name"),
                                          [l[1]]))
            else:
                print "WARNING unexpected line in jobs section:",repr(line)
            
mainConv.add_conversion(Conv("JOBS",None,special_convert=jobs_conv))

def ot_conv(oldInput,oldSect,newInput,new_sect,conv):
    if not new_sect.args:
        new_sect.args=["OFF"]
mainConv.add_conversion(Conv("OTSCF","FORCE_EVAL%DFT%SCF%OT",
                        special_convert=ot_conv))

def top_conv(oldInput,oldSect,newInput,new_sect,conv):
    if new_sect.keywords.has_key("forcefield"):
        del new_sect.keywords['forcefield']
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("TOPOL","FORCE_EVAL%SUBSYS%TOPOLOGY",
                             special_convert=top_conv))
mainConv.add_conversion(Conv("TOPOLOGY","FORCE_EVAL%SUBSYS%TOPOLOGY",
                             special_convert=top_conv))

def mc_conv(oldInput,oldSect,newInput,new_sect,conv):
    for k in new_sect.keywords_list:
        if k.name[-1:]==':':
            k.name=k.name[:-1]
        if len(k.values)>1:
            if k.values[1][0]=="(":
                k.values=[k.values[0]]
mainConv.add_conversion(Conv("MC","MOTION%MC",
                             special_convert=mc_conv))

def glob_conv(oldInput,oldSect,newInput,new_sect,conv):
    if (new_sect.keywords.has_key("iolevel") and
        new_sect.keywords.has_key("print_level")):
        del new_sect.keywords["iolevel"]
        new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("GLOBAL","GLOBAL",special_convert=glob_conv))

mainConv.add_conversion(Conv("GEOOPT","MOTION%GEOOPT"))

mainConv.add_conversion(Conv("MOTION","MOTION"))

def ew_conv(oldInput,oldSect,newInput,new_sect,conv):
    k=oldSect.keywords['ewald_type'][0]
    ew_type=k.values[0].lower()
    new_sect.add_keyword(Keyword(conv.upcase(k.name),[k.values[0]]))
    ew_params={'spme':{0:'alpha',1:'gmax',2:'o_spline'},
              'pme':{0:'alpha',1:'ns_max',2:'epsilon'},
              'ewald':{0:'alpha',1:'gmax'}}
    if oldSect.keywords.has_key('ewald_param'):
        k=oldSect.keywords['ewald_param'][0]
        ew_param=ew_params[ew_type]
        for idx in ew_param.keys():
            if len(k.values)>idx:
                new_sect.add_keyword(Keyword(conv.upcase(ew_param[idx]),
                                             [k.values[idx]]))
mainConv.add_conversion(Conv("EWALD","FORCE_EVAL%MM%POISSON_MM%EWALD",
                            merge=0,special_convert=ew_conv))

def coord_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    new_sect.pre_comments=oldSect.pre_comments
    new_sect.post_comments=oldSect.post_comments
    new_sect.raw_lines=map(lambda x: "   "+x,oldSect.raw_lines)
mainConv.add_conversion(Conv("COORD","FORCE_EVAL%SUBSYS%COORD",merge=0,
                             special_convert=coord_conv))
mainConv.add_conversion(Conv("COORDS","FORCE_EVAL%SUBSYS%COORD",merge=0,
                             special_convert=coord_conv))

def ff_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    new_sect.pre_comments=oldSect.pre_comments
    new_sect.post_comments=oldSect.post_comments
    l_nr=0
    nl=len(oldSect.raw_lines)
    while 1:
        if l_nr>=nl: break
        line=oldSect.raw_lines[l_nr]
        ll=line.strip().lower()
        if ll=="charges":
            while 1:
                l_nr+=1
                if l_nr>=nl: break
                line=oldSect.raw_lines[l_nr]
                if line.split()[0].lower()=="end":break
                ch=Section(conv.upcase("charge"))
                ch.add_keyword(Keyword(conv.upcase("atom"),
                                       values=[line.split()[0]]))
                ch.add_keyword(Keyword(conv.upcase("charge"),
                                       values=[line.split()[1]]))
                new_sect.add_subsection(ch)
        elif ll=="bends":
            while 1:
                l_nr+=1
                if l_nr>=nl: break
                line=oldSect.raw_lines[l_nr]
                if line.split()[0].lower()=="end":break
                ch=Section(conv.upcase("bend"))
                ch.add_keyword(Keyword(conv.upcase("atoms"),values=line.split()[1:4]))
                ch.add_keyword(Keyword(conv.upcase("k"),values=[line.split()[4]]))
                ch.add_keyword(Keyword(conv.upcase("theta0"),values=[line.split()[5]]))
                new_sect.add_subsection(ch)
        elif ll=="bonds":
            while 1:
                l_nr+=1
                if l_nr>=nl: break
                line=oldSect.raw_lines[l_nr]
                ll=line.split()
                if ll[0].lower()=="end":break
                ch=Section(conv.upcase("bond"))
                ch.add_keyword(Keyword(conv.upcase("atoms"),values=ll[1:3]))
                if ll[0].lower()=="harmonic":
                    ch.add_keyword(Keyword(conv.upcase("k"),values=[ll[3]]))
                    ch.add_keyword(Keyword(conv.upcase("r0"),values=[line.split()[4]]))
                elif ll[0].lower()=="quartic":
                    ch.add_keyword(Keyword(conv.upcase("k"),values=ll[3:6]))
                    ch.add_keyword(Keyword(conv.upcase("r0"),values=[line.split()[6]]))
                else:
                    print "WARNING unknown bond type in forcefield section:",ll[0]
                new_sect.add_subsection(ch)
        elif ll.split()[0]=="parmfile":
            new_sect.add_keyword(Keyword("parmfile",[line.split()[2]]))
            new_sect.add_keyword(Keyword("parmtype",[line.split()[1]]))
        elif ll.split()[0]=="ei_scale":
            new_sect.add_keyword(Keyword("ei_scale14",[line.split()[1]]))
        elif ll.split()[0]=="vdw_scale":
            new_sect.add_keyword(Keyword("vdw_scale14",[line.split()[1]]))
        elif ll.split()[0]=="rcut_nb":
            new_sect.add_keyword(Keyword("rcut_nb",[line.split()[1]]))
        elif ll.split()[0]=="nonbonded":
            ss=Section(conv.upcase("nonbonded"))
            new_sect.add_subsection(ss)
            f_data={'lennard-jones':{3:'epsilon',4:'sigma',5:'rcut'},
                    'bmhft':{3:'rcut'},'eam':{3:'parmfile'},'ipbv':{3:'rcut'},
                    'williams':{3:'a',4:'b',5:'c',6:'rcut'}}
            while 1:
                l_nr+=1
                if l_nr>=nl: break
                line=oldSect.raw_lines[l_nr]
                if line.split()[0].lower()=="end":break
                l=line.split()
                sname=l[0].lower()
                ch=Section(conv.upcase(l[0]))
                ch.add_keyword(Keyword("atoms",l[1:3]))
                for idx in f_data[sname]:
                    kname=f_data[sname][idx]
                    ch.add_keyword(Keyword(conv.upcase(kname),values=[l[idx]]))
                ss.add_subsection(ch)
        else:
            print "WARNING ignoring line ",repr(line),"in forcefield section"
        l_nr+=1          
mainConv.add_conversion(Conv("FORCE_FIELD","FORCE_EVAL%MM%FORCEFIELD",merge=0,
                             special_convert=ff_conv))

def io_conv(oldInput,oldSect,newInput,new_sect,conv):
    for k in oldSect.keywords_list:
        kname=k.name.lower()
        if (kname=='basis_set_file' or
            kname=='potential_file' or kname=='restart_file'):
            s=guaranteePath(newInput,conv.upcase("FORCE_EVAL%DFT"))
            s[-1].add_keyword(Keyword(conv.upcase(k.name+"_name"),
                                      [k.values[0]]))
        elif (kname=='mm_potential_file'):
            s=guaranteePath(newInput,conv.upcase("FORCE_EVAL%QMMM"))
            s[-1].add_keyword(Keyword(conv.upcase(k.name+"_name"),
                                      [k.values[0]]))
        else:
            print "WARNING dropping &IO keyword ",k.name
mainConv.add_conversion(Conv("IO",None,special_convert=io_conv))

if __name__=="__main__":
    if len(sys.argv)>3 or len(sys.argv)<2:
        print os.path.basename(sys.argv[0])+" input_file [output_file]"
        sys.exit(1)
    inF=file(sys.argv[1])
    parser=Parser(inF)
    outF=sys.stdout
    if len(sys.argv)==3:
        outF=file(sys.argv[2],"w")
    oldInput=InputFile(inF.name)
    oldInput.parse_file(parser)
    newInput=InputFile(outF.name)
    mainConv.convert(oldInput,newInput)
    newInput.write(outF)
