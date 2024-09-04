<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output indent="no" method="text" omit-xml-declaration="yes" version="2.0"/>
<xsl:template match="/">
" Vim syntax file
" Language: CP2K input files
" Maintainer: based on work by Alin M Elena available at https://gitorious.org/cp2k-vim
" Revision: 13th of September 2010
" History:
" - XSLT dump and improved syntax highlighting (10.12.2013, Matthias Krack)
" - Folding and automatic indentation added (13.12.2013, Matthias Krack)
" - Remove folding since it overrides user's defaults (18.11.2016, Patrick Seewald)
" CP2K-Version: <xsl:value-of select="/CP2K_INPUT/CP2K_VERSION"/> (<xsl:value-of select="/CP2K_INPUT/COMPILE_REVISION"/>)

if exists("b:current_syntax")
   finish
endif
let b:current_syntax = "cp2k"

syn case ignore

"----------------------------------------------------------------/
" CP2K numbers used in blocks
"----------------------------------------------------------------/
" Regular int like number with - + or nothing in front
" e.g. 918, or -49

syn match cp2kNumber '\d\+'
syn match cp2kNumber '[-+]\d\+'

" Floating point number with decimal no E or e (+,-)
" e.g. 0.198781, or -3.141592

syn match cp2kNumber '\d\+\.\d*'
syn match cp2kNumber '[-+]\d\+\.\d*'

" Floating point like number with E and no decimal point (+,-)
" e.g. 3E+9, 3e+09, or -3e+02

syn match cp2kNumber '[-+]\=\d[[:digit:]]*[eE][\-+]\=\d\+'
syn match cp2kNumber '\d[[:digit:]]*[eE][\-+]\=\d\+'

" Floating point like number with E and decimal point (+,-)
" e.g. -3.9188e+09, or 0.9188E-93

syn match cp2kNumber '[-+]\=\d[[:digit:]]*\.\d*[eE][\-+]\=\d\+'
syn match cp2kNumber '\d[[:digit:]]*\.\d*[eE][\-+]\=\d\+'

"----------------------------------------------------------------/
" CP2K Boolean entities
"----------------------------------------------------------------/

syn keyword cp2kBool true false T F yes no

"-----------------------------------------------------------------/
" CP2K comments
"-----------------------------------------------------------------/

syn keyword cp2kTodo TODO FIXME NOTE REMARK
syn match cp2kComment "#.*$" contains=cp2kTodo
syn match cp2kComment "!.*$" contains=cp2kTodo

"----------------------------------------------------------------/
" CP2K predefined constants
"----------------------------------------------------------------/
<xsl:for-each-group select="//ITEM" group-by="NAME">
<xsl:sort select="NAME"/>
syn keyword cp2kConstant <xsl:value-of select="NAME"/>
</xsl:for-each-group>

"----------------------------------------------------------------/
" CP2K sections
"----------------------------------------------------------------/
<xsl:for-each-group select="//SECTION" group-by="NAME">
<xsl:sort select="NAME"/>
syn keyword cp2kSection <xsl:value-of select="NAME"/>
</xsl:for-each-group>
syn keyword cp2kSection END
syn keyword cp2kSection ENDIF
syn match cp2kBegSection '^\s*&amp;\w\+' contains=cp2kSection
syn match cp2kEndSection '^\s*&amp;END\s*\w\+' contains=cp2kSection

"----------------------------------------------------------------/
" CP2K default keyword names
"----------------------------------------------------------------/
<xsl:for-each-group select="//KEYWORD" group-by="NAME[@type='default']">
<xsl:sort select="NAME[@type='default']"/>
<xsl:if test="NAME[@type='default']!=''">
syn keyword cp2kKeyword <xsl:value-of select="NAME[@type='default']"/>
</xsl:if>
</xsl:for-each-group>

"----------------------------------------------------------------/
" CP2K alias keyword names
"----------------------------------------------------------------/
<xsl:for-each-group select="//KEYWORD" group-by="NAME[@type='alias'][1]">
<xsl:sort select="NAME[@type='alias'][1]"/>
<xsl:if test="NAME[@type='alias']!=''">
syn keyword cp2kKeyword <xsl:value-of select="NAME[@type='alias']"/>
</xsl:if>
</xsl:for-each-group>

"-----------------------------------------------------------------/
" CP2K preprocessing directives
"-----------------------------------------------------------------/

syn keyword cp2kPreProc ENDIF FFTYPE IF INCLUDE PRINT SET XCTYPE

"-----------------------------------------------------------------/
" CP2K strings
"-----------------------------------------------------------------/

syn region cp2kString matchgroup=cp2kStringDelimiter start=+"+ end=+"+
syn region cp2kString matchgroup=cp2kStringDelimiter start=+'+ end=+'+
syn region cp2kString matchgroup=cp2kStringDelimiter start=+`+ end=+`+

"----------------------------------------------------------------------------/
" Final setup
"----------------------------------------------------------------------------/

setlocal autoindent
setlocal expandtab
setlocal iskeyword+=-

"----------------------------------------------------------------------------/
" Indentation support for CP2K input syntax
"----------------------------------------------------------------------------/

if exists("b:did_indent")
   finish
endif
let b:did_indent = 1

setlocal indentexpr=GetCp2kIndent()
setlocal indentkeys+=0=~&amp;END,0=~@ENDIF,0=#,0=@
setlocal nosmartindent

if exists("*GetCp2kIndent")
   finish
endif

function! GetCp2kIndent()
   let lnum = prevnonblank(v:lnum - 1)
   if lnum == 0
      return 0
   endif
   let ind = indent(lnum)
   let line = getline(lnum)
   if line =~ '^\s*\c\%(&amp;\|@IF\)'
      let ind += &amp;shiftwidth
      if line =~ '^\s*\c\%(&amp;END\|@ENDIF\)\>'
         let ind -= &amp;shiftwidth
      endif
   endif
   let line = getline(v:lnum)
   if line =~ '^\s*\c\%(&amp;END\|@ENDIF\)\>'
      let ind -= &amp;shiftwidth
   endif
   return ind
endfunction

"----------------------------------------------------------------------------/
" CP2K keyword highlighting rules
"----------------------------------------------------------------------------/

hi def link cp2kComment         Comment
hi def link cp2kConstant        Constant
hi def link cp2kTodo            Todo
hi def link cp2kBool            Boolean
hi def link cp2kNumber          Number
hi def link cp2kKeyword         Keyword
hi def link cp2kSection         Structure
hi def link cp2kPreProc         PreProc
hi def link cp2kString          String
hi def link cp2kStringDelimiter Delimiter
</xsl:template>
</xsl:stylesheet>
