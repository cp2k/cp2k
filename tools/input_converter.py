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
                            print "        having keywords, probably something is wrong"
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
    if new_sect.keywords.has_key('functional'):
        p=guaranteePath(newInput,conv.upcase("force_eval%dft%xc%xc_functional"))
        del p[-1].auto_g
        p[-1].args=new_sect.keywords['functional'][0].values
        del new_sect.keywords['functional']
        new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("DFT","FORCE_EVAL%DFT",special_convert=dft_conv))

mainConv.add_conversion(Conv("FORCE_EVAL","FORCE_EVAL"))

def cp2k_conv(oldInput,oldSect,newInput,new_sect,conv):
    if new_sect.keywords.has_key('iolevel'):
        del new_sect.keywords['iolevel']
        new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("CP2K","GLOBAL",special_convert=cp2k_conv))

def md_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    for k in oldSect.keywords_list:
        easyT={'ensem':'ensemble','steps':'steps','tempe':'temperature',
               'resta':'restart','temp_':'temp_tol','times':'timestep'}
        noseT={'lengt':'length','yoshi':'yoshida','timec':'timecon',
               'mts:':'mts'}
        if easyT.has_key(k.name[:5].lower()):
            new_sect.add_keyword(Keyword(conv.upcase(easyT[k.name[:5].lower()]),
                                         values=[conv.upcase(k.values[0])]))
        elif noseT.has_key(k.name[:5].lower()):
            s=guaranteePath(newInput,conv.upcase("motion%md_new%nose"))
            if hasattr(s[-1],'auto_g'):
                del s[-1].auto_g
            s[-1].add_keyword(Keyword(conv.upcase(noseT[k.name[:5].lower()]),
                                      values=[conv.upcase(k.values[0])]))
        elif k.name[:5].lower()=='shake':
            s=guaranteePath(newInput,conv.upcase("subsys%constraint"))
            s[-1].add_keyword(Keyword(conv.upcase("shake"),
                                      values=[conv.upcase(k.values[0])]))
        elif k.name[:5].lower()!='const' and k.name[:5].lower()!='nose_':
            print "WARNING ignoring unknown keyword",k.name,"in section MD"
mainConv.add_conversion(Conv("MD","MOTION%MD_NEW",special_convert=md_conv,merge=0))

def qs_conv(oldInput,oldSect,newInput,new_sect,conv):
    s=guaranteePath(newInput,conv.upcase("FORCE_EVAL%DFT%MGRID"))
    toMg=['cutoff','ngrids','progression_factor',
          'commensurate','realspace','rel_cutoff']
    for k in new_sect.keywords_list:
        if k.name.lower() in toMg:
            s[-1].add_keyword(k)
            del new_sect.keywords[k.name.lower()]
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("QS","FORCE_EVAL%DFT%QS",special_convert=qs_conv))

def scf_conv(oldInput,oldSect,newInput,new_sect,conv):
    for k in new_sect.keywords_list:
        if k.name.lower()=="guess":
            k.name=conv.upcase("scf_guess")
            new_sect.keywords[k.name.lower()]=[k]
            del new_sect.keywords["guess"]
        if k.name.lower()=="ot":
            s=guaranteePath(new_sect,conv.upcase("ot"))
            s[-1].args=["ON"]
            del new_sect.keywords["ot"]
    new_sect.rebuild_keyword_list()
mainConv.add_conversion(Conv("SCF","FORCE_EVAL%DFT%SCF",
                             special_convert=scf_conv))

def print_conv(oldInput,oldSect,newInput,new_sect,conv):
    if oldSect.args:
        s=guaranteePath(newInput,conv.upcase("global"))
        s[-1].add_keyword(Keyword(conv.upcase("print_level"),
                                  [conv.upcase(oldSect.args[0])]))
    for k in oldSect.keywords_list:
        print "WARNING dropping printkey",k.name
mainConv.add_conversion(Conv("PRINT",None,special_convert=print_conv))

mainConv.add_conversion(Conv("KIND","SUBSYS%KIND"))

mainConv.add_conversion(Conv("CELL","SUBSYS%CELL"))

def coord_conv(oldInput,oldSect,newInput,new_sect,conv):
    del new_sect.auto_g
    new_sect.pre_comments=oldSect.pre_comments
    new_sect.post_comments=oldSect.post_comments
    new_sect.raw_lines=map(lambda x: "   "+x,oldSect.raw_lines)
mainConv.add_conversion(Conv("COORD","SUBSYS%COORD",merge=0,
                             special_convert=coord_conv))

if __name__=="__main__":
    if len(sys.argv)>3 or len(sys.argv)<2:
        print os.path.basename(sys.argv[0])+" input_file [output_file]"
        sys.exit(1)
    inF=file(sys.argv[1])
    parser=Parser(inF)
    outF=sys.stdout
    if len(sys.argv)==3:
        outF=file(sys.argv[2])
    oldInput=InputFile(inF.name)
    oldInput.parse_file(parser)
    newInput=InputFile(outF.name)
    mainConv.convert(oldInput,newInput)
    newInput.write(outF)
