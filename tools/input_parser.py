#!env python
import sys, os, os.path, StringIO

class Parser:
    "represent a file to parse"
    def __init__(self,f,lines=[],lineN=0):
        self.f=f
        self.lines=lines
        self.lineN=lineN
    def readline(self):
        self.lineN+=1
        if self.lines:
            return self.lines.pop()
        line=self.f.readline()
        if not line: self.lineN-=1
        return line
    def readInputLine(self):
        coreLines=[]
        lines=""
        comment=""
        while 1:
            line=self.readline()
            if not line: break
            lines+=line
            commStart=line.find('#')
            if commStart>=0:
                comment+=" "+line[commStart+1:].strip(" \t\n")
                coreLines.append(line[:commStart].strip())
            else:
                coreLines.append(line.strip(" \t\n"))
            if not coreLines[-1] or coreLines[-1][-1]!="\\": break
            coreLines[-1]=coreLines[-1][-1].strip(" \t\n")
        return (" ".join(coreLines),comment,lines)
    def pushback(self,lines):
        if lines:
            for line in lines.splitlines(1):
                self.lineN-=1
                self.lines.append(line)
    def file_pos(self):
        line=self.readline()
        self.pushback(line)
        return ("%s at line %i: %s"%(repr(self.f.name),self.lineN,repr(line)))
    def old_pos(self,line=None):
        return ("%s at line %i: %s"%(repr(self.f.name),self.lineN-1,repr(line)))
    def close(self):
        self.f.close()

class Keyword:
    def __init__(self,name,values,comment=None,pre_comments=None):
        self.name=name
        if values:
            self.values=values
        else:
            self.values=[]
        self.comment=comment
        if pre_comments:
            self.pre_comments=pre_comments
        else:
            self.pre_comments=[]
    def write(self,outF,indentLevel=0):
        for comment in self.pre_comments:
            outF.write("  "*(indentLevel-1)+"#"+comment+"\n")
        outF.write("  "*indentLevel+self.name)
        ll=2*indentLevel+len(self.name)
        for val in self.values:
            if ll+len(val)+1>80:
                outF.write("\\\n "+("  "*indentLevel))
                ll=2*indentLevel+1
            outF.write(" "+val)
            ll+=len(val)+1
        if self.comment:
            outF.write(" #"+self.comment)
        outF.write("\n")
    def __str__(self):
        outS=StringIO.StringIO()
        self.write(outS)
        return outS.getvalue()
    
class Section:
    def __init__(self,name,args=None,subsections=None,raw_lines=None,
                 keywords_list=None,pre_comments=None,post_comments=None):
        self.name=name
        if args:
            self.args=args
        else:
            self.args=[]
        if subsections:
            self.subsections=subsections
        else:
            self.subsections={}
        if raw_lines:
            self.raw_lines=raw_lines
        else:
            self.raw_lines=[]
        if keywords_list:
            self.keywords_list=keywords_list
        else:
            self.keywords_list=[]
        self.keywords={}
        for kw in self.keywords_list:
            if (not self.keywords.has_key(kw.name.lower())):
                self.keywords[kw.name.lower()]=[]
            self.keywords[kw.name.lower()].append(kw)
        if pre_comments:
            self.pre_comments=pre_comments
        else:
            self.pre_comments=[]
        if post_comments:
            self.post_comments=post_comments
        else:
            self.post_comments=[]
    def add_input_line(self,line,comment,lines):
        self.raw_lines.append(lines)
        if line:
            l=line.split()
            kw=Keyword(name=l[0],values=l[1:],comment=comment,
                       pre_comments=self.post_comments)
            self.post_comments=[]
            self.add_keyword(kw)
        elif comment:
            self.post_comments.append(comment)
    def add_keyword(self,kw):
        self.keywords_list.append(kw)
        if (not self.keywords.has_key(kw.name.lower())):
            self.keywords[kw.name.lower()]=[]
        self.keywords[kw.name.lower()].append(kw)
    def add_subsection(self,subS):
        if not self.subsections.has_key(subS.name.lower()):
            self.subsections[subS.name.lower()]=[]
        self.subsections[subS.name.lower()].append(subS)
    def write(self,outF,indentLevel=0):
        indentC="  "*(indentLevel-1)+"#"
        for comment in self.pre_comments:
            outF.write(indentC)
            outF.write(comment)
            outF.write("\n")
        if self.name==None: return
        indentS="  "*indentLevel
        outF.write(indentS)
        outF.write("&"+self.name)
        for arg in self.args:
            outF.write(" "+arg)
        outF.write("\n")
        if self.raw_lines:
            outF.writelines(self.raw_lines)
        else:
            for kw in self.keywords_list:
                kw.write(outF,indentLevel+1)
        sNames=self.subsections.keys()
        sNames.sort()
        for sName in sNames:
            for subS in self.subsections[sName]:
                subS.write(outF,indentLevel=indentLevel+1)
        indentC=indentS+"#"
        for comment in self.post_comments:
            outF.write(indentC)
            outF.write(comment)
            outF.write("\n")
        outF.write(indentS)
        if (self.name.lower()==self.name):
            outF.write("&end ")
        else:
            outF.write("&END ")
        outF.write(self.name)
        outF.write("\n")
    def remove_raw_lines(self,recursive=1):
        self.raw_lines=[]
        if (recursive):
            for subSs in self.subsections.values():
                for subS in subSs:
                    subS.remove_raw_lines(recursive)
    def __str__(self):
        outS=StringIO.StringIO()
        self.write(outS)
        return outS.getvalue()
    def sort_keywords(self,recursive=0):
        self.keywords_list.sort(lambda x,y:cmp(x.name.lower(),y.name.lower()))
        if (recursive):
            for subSs in self.subsections.values():
                for subS in subSs:
                    subS.keywords_list.sort(lambda x,y:cmp(x[0].name.lower(),y[0].name.lower()))
    def rebuild_keyword_list(self,recursive=0):
        self.keywords_list=[]
        for kwrds in self.keywords.itervalues():
            for kwrd in kwrds:
                self.keywords_list.append(kwrd)
        self.sort_keywords(recursive=0)
        if recursive:
            for subSs in self.subsections.itervalues():
                for subS in subSs:
                    subS.rebuild_keyword_list(self,recursive=recursive)
        
def read_section(f):
    while 1:
        (coreLine,comments,lines)=f.readInputLine()
        if not lines or not lines.isspace(): break
    preSectionComments=[]
    while lines and (not coreLine or coreLine[0] != "&"):
        if comments:
            preSectionComments.append(comments)
        if coreLine:
            print "ignoring out of section line ",f.old_pos(coreLine)
        (coreLine,comments,lines)=f.readInputLine()
    if not lines:
        if preSectionComments:
            return Section(None,pre_comments=preSectionComments)
        else:
            return None
    l=coreLine.split()
    sect_name=l[0][1:]
    sect_args=l[1:]
    if sect_name.lower()=="end": raise SyntaxError("found dangling &END in %s"%f.old_pos(line))
    sect=Section(sect_name,sect_args,pre_comments=preSectionComments)
    while 1:
        (coreLine,comments,lines)=f.readInputLine()
        if not lines:
            raise SyntaxError("unterminated section %s in %s"%
                              (sect_name,f.old_pos(line)))
        if coreLine and coreLine[0]=='&':
            l=coreLine.split()
            if l[0][1:].lower()=="end":
                if len(l)>1 and l[-1].lower()!=sect_name.lower():
                    raise SyntaxError("end section mismatch %s vs %s at line %i"%
                                      (sect_name,l[-1],f.lineN))
                else:
                    return sect
            else:
                f.pushback(lines)
                sect.add_subsection(read_section(f))
        else:
            sect.add_input_line(coreLine,comments,lines)

class InputFile:
    def __init__(self,filename,pre_comments=None,post_comments=None,subsections=None):
        self.filename=filename
        if pre_comments:
            self.pre_comments=pre_comments
        else:
            self.pre_comments=[]
        if post_comments:
            self.post_comments=post_comments
        else:
            self.post_comments=[]
        if subsections:
            self.subsections=subsections
        else:
            self.subsections={}
    def parse_file(self,f):
        while 1:
            (coreLine,comments,lines)=f.readInputLine()
            if not lines: break
            if coreLine or lines.isspace():
                f.pushback(lines)
                break
            self.pre_comments.append(comments)
        while 1:
            sect=read_section(f)
            if not sect: break
            if sect.name==None:
                self.post_comments+=sect.pre_comments
                break
            self.add_subsection(sect)
    def add_subsection(self,sect):
        if not self.subsections.has_key(sect.name.lower()):
            self.subsections[sect.name.lower()]=[]
        self.subsections[sect.name.lower()].append(sect)
    def write(self,outF):
        for comment in self.pre_comments:
            outF.write("#")
            outF.write(comment)
            outF.write("\n")
        outF.write("\n")
        sNames=self.subsections.keys()
        sNames.sort()
        for sName in sNames:
            for sect in self.subsections[sName]:
                sect.write(outF)
        for comment in self.post_comments:
            outF.write("#")
            outF.write(comment)
            outF.write("\n")
    def remove_raw_lines(self):
        for ss in self.subsections.values():
            for s in ss:
                s.remove_raw_lines(1)
        
if __name__=="__main__":
    if len(sys.argv)>3 or len(sys.argv)<2:
        print os.path.basename(sys.argv[0])+" input_file [output_file]"
        sys.exit(1)
    inF=file(sys.argv[1])
    parser=Parser(inF)
    outF=sys.stdout
    if len(sys.argv)==3:
        outF=file(sys.argv[2])
    parsedFile=InputFile(inF.name)
    parsedFile.parse_file(parser)
    print "*** parsing complete ***"
    parsedFile.remove_raw_lines()
    parsedFile.write(outF)
    print "*** done ***"
