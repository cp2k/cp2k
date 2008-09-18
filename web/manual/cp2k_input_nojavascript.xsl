<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">

<xsl:output method="html" indent="yes" name="html"/>

<xsl:template match="/CP2K_INPUT">
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
   <xsl:call-template name="head"/>
   <title>CP2K input reference</title>
  </head>
  <body>
   <h1 align="center">CP2K input reference</h1>
   <h2>Version information</h2>
   <p>
    This HTML manual was generated automatically from a CP2K executable
    compiled on <xsl:value-of select="COMPILE_DATE"/> using the
    --xml command line option. Thus the manual describes exactly this
    version of the code. The latest CVS log file entry found was
    <xsl:value-of select="COMPILE_LASTCVS"/>.
   </p>
   <p>
    <xsl:call-template name="searchform"/>
   </p>
   <h2>Journal papers</h2>
   <p>
    <a href="references.html">List of references</a> cited in the CP2K input manual.
   </p>
   <h2>Internal input preprocessor</h2>
   <p>
    Before the input is parsed, the input is run through a simple internal preprocessor.
    The preprocessor recognizes the following directives independent of capitalization:
   </p>
   <dl>
    <dt><b>@INCLUDE 'filename.inc'</b></dt>
    <dd>
     The file referenced by <i>filename.inc</i> is included into the input file and parsed.
     Recursive inclusions are not allowed and the files have to exist in the current working
     directory. There can be only one @INCLUDE statement per line. Single or double quotes
     have to be used if the filename contains blanks.
    </dd>
    <dt><b>@SET VAR value</b></dt>
    <dd>
     Assigns the text <i>value</i> to the preprocessing variable <i>VAR</i>. <i>value</i>
     is the text following <i>VAR</i> with the outer whitespace removed. The variable can be
     recalled with a <i>${VAR}</i> statement. There can be only one @SET statement per line.
    </dd>
    <dt><b>${VAR}</b></dt>
    <dd>
     Expand the variable <i>VAR</i>. The text <i>${VAR}</i> is replaced with the value assigned
     to <i>VAR</i> in the last @SET directive. There can be multiple variable statements
     per line. The expansion process is repeated until no more variables are found.
    </dd>
     <dt><b>@IF / @ENDIF</b></dt>
     <dd>
      Conditional block. The text from the @IF line up to the next line with a valid "
      @ENDIF is skipped, if the expression following @IF resolves to <i>false</i>.
      Available expressions are lexical comparisons for equality '==' or inequality '/='.
      If none of the two operators are found, a '0' or whitespace resolves to <i>false</i>
      while any text resolves to <i>true</i>. @IF/@ENDIF blocks cannot be nested and
      cannot span across files. There can be only one test (== or /=) per @IF statement.
    </dd>
   </dl>
   <h2>Input structure</h2>
    All sections and keywords that can be part of a CP2K input file are shown
    with their allowed nestings. A detailed description can be obtained by
    clicking on the section links. The links in the detailed descriptions switch
    back to the corresponding index entries. In this way a toggling between the
    index and the detailed description is feasible.
   <h2>Index of all input sections</h2>
   <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
   <ul class="disc">
    <li>
     <a href="CP2K_INPUT.html" id="CP2K_INPUT.html">CP2K_INPUT</a>
    </li>
    <xsl:call-template name="describe_sections">
     <xsl:with-param name="path" select="'CP2K_INPUT'"/>
     <xsl:with-param name="root" select="'../'"/>
    </xsl:call-template>
   </ul>
   <p>
    <a href="index_javascript.html">Collapsible section tree</a> (CAUTION: very long loading time, ignore browser warnings, just continue)
   </p>
   <xsl:call-template name="footer">
    <xsl:with-param name="root" select="''"/>
   </xsl:call-template>
  </body>
  <xsl:result-document href="CP2K_INPUT.html" method="html" indent="yes" format="html">
   <html>
    <head>
     <xsl:call-template name="head"/>
     <title>The CP2K project: CP2K input file</title>
    </head>
    <body>
     <xsl:call-template name="header">
      <xsl:with-param name="root" select="''"/>
     </xsl:call-template>
     <h2><a href="index.html#CP2K_INPUT.html">CP2K input file</a></h2>
     <ul class="none">
      <li>
       <em>Input file of CP2K</em>
      </li>
     </ul>
     <ul class="none">
      <li>
       Section path:
       <big class="uctt">
        <a href="CP2K_INPUT.html">CP2K_INPUT</a>
       </big>
      </li>
     </ul>
     <ul class="none">
      <li>
       This section cannot be repeated and can be optional.
      </li>
     </ul>
     <ul class="none">
      <h3>Subsections</h3>
      <xsl:for-each select="SECTION">
       <xsl:sort select="NAME"/>
       <ul class="disc">
        <li>
         <a href="{concat('CP2K_INPUT/',string(NAME),'.html')}"><xsl:value-of select="NAME"/></a>
        </li>
       </ul>
      </xsl:for-each>
     </ul>
     <ul class="none">
      <h3>Keywords</h3>
      <ul class="disc">
       <li>
        none
       </li>
      </ul>
     </ul>
     <xsl:call-template name="footer">
      <xsl:with-param name="root" select="''"/>
     </xsl:call-template>
    </body>
   </html>
  </xsl:result-document>
 </html>
</xsl:template>

<xsl:template name="head">
 <meta name="description" content="CP2K"/>
 <meta name="keywords" contents="scientific,computing,chemistry,physics,documentation,help,manual,Fortran,parallel"/>
 <style type="text/css">
  body {background-color: #ffffff}
  big.tt {font-family: monospace; font-size: 100%}
  big.uctt {font-family: monospace; font-size: 100%; text-transform: uppercase}
  p.uctt {font-family: monospace; text-transform: uppercase}
  table.default {table-layout: fixed; width: 100%}
  td.l {width: 25%}
  td.r {width: 75%}
  ul.circle {list-style-type: circle}
  ul.disc {list-style-type: disc}
  ul.none {list-style-type: none}
  ul.square {list-style-type: square}
 </style>
 <script src="collapsibleList.js" type="text/javascript" language="javascript1.2"/>
</xsl:template>

<xsl:template name="header">
 <xsl:param name="root"/>
 <table class="default">
  <tr>
   <td align="left">
    Back to the <a href="{$root}index.html">main page</a> of this manual
   </td>
   <td align="right">
    <xsl:call-template name="searchform"/>
   </td>
  </tr>
 </table>
 <hr/>
</xsl:template>

<xsl:template name="footer">
 <xsl:param name="root"/>
 <hr/>
 <table class="default">
  <tr>
   <td align="left">
    Back to the <a href="{$root}index.html">main page</a> of this manual or the <a href="http://cp2k.berlios.de">CP2K home page</a>
   </td>
   <td align="right">
    (Last update:
    <xsl:value-of select="day-from-dateTime(current-dateTime())"/>.<xsl:value-of select="month-from-dateTime(current-dateTime())"/>.<xsl:value-of select="year-from-dateTime(current-dateTime())"/>)
   </td>
  </tr>
 </table>
</xsl:template>

<xsl:template name="searchform">
 <form method="get" action="http://www.google.com/search">
  <input type="text" name="q" maxlength="255"/>
  <input type="submit" value="Search this manual (Google)"/>
  <input type="hidden" name="domains" value="http://cp2k.berlios.de/manual/"/>
  <input type="radio" style="visibility:hidden" name="sitesearch" value="http://cp2k.berlios.de/manual/" checked="checked"/>
 </form>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:param name="path"/>
 <xsl:param name="root"/>
 <ul class="disc">
  <xsl:for-each select="SECTION">
   <xsl:sort select="NAME"/>
   <xsl:variable name="section" select="string(NAME)"/>
   <xsl:variable name="local_path" select="concat($path,'/',string(NAME))"/>
   <xsl:variable name="section_filename" select="concat($local_path,'.html')"/>
   <li>
    <a href="{$section_filename}" id="{$section_filename}"><xsl:value-of select="NAME"/></a>
   </li>
   <xsl:result-document href="{$section_filename}" method="html" indent="yes" format="html">
   <html>
    <head>
     <xsl:call-template name="head"/>
     <title>The CP2K project: input section <xsl:value-of select="NAME"/></title>
    </head>
    <body>
     <xsl:call-template name="header">
      <xsl:with-param name="root" select="$root"/>
     </xsl:call-template>
     <h2><a href="{$root}index.html#{$section_filename}">Section <xsl:value-of select="NAME"/></a></h2>
     <xsl:if test="string-length(DESCRIPTION) > 0">
      <ul class="none">
       <li>
        <em>
         <xsl:value-of select="DESCRIPTION"/>
        </em>
       </li>
      </ul>
     </xsl:if>
     <ul class="none">
      <li>
       Section path:
       <big class="uctt">
        <xsl:call-template name="link_section_path">
         <xsl:with-param name="string" select="$local_path"/>
         <xsl:with-param name="separator" select="'/'"/>
         <xsl:with-param name="root" select="$root"/>
        </xsl:call-template>
       </big>
      </li>
     </ul>
     <ul class="none">
      <li>
       This section can<xsl:if test="@repeats = 'no'">not</xsl:if> be repeated
       and can<xsl:if test="@required = 'yes'">not</xsl:if> be optional.
      </li>
     </ul>
     <xsl:if test="count(REFERENCE) > 0">
      <ul class="none">
       <li>
        This section cites the following reference<xsl:if test="count(REFERENCE) > 1">s</xsl:if>:
        <xsl:for-each select="REFERENCE">
         <xsl:sort select="NAME"/>
         [<a href="{$root}references.html#reference_{string(NUMBER)}"><xsl:value-of select="NAME"/></a>]
        </xsl:for-each>
       </li>
      </ul>
     </xsl:if>
     <ul class="none">
      <li>
       <h3>Subsections</h3>
       <xsl:choose>
        <xsl:when test="count(SECTION) > 0">
         <ul class="disc">
          <xsl:for-each select="SECTION">
           <xsl:sort select="NAME"/>
           <xsl:variable name="subsection_filename" select="concat($section,'/',string(NAME),'.html')"/>
           <li>
            <a href="{$subsection_filename}"><xsl:value-of select="NAME"/></a>
           </li>
          </xsl:for-each>
         </ul>
        </xsl:when>
        <xsl:otherwise>
         <ul class="disc">
          <li>
           none
          </li>
         </ul>
        </xsl:otherwise>
       </xsl:choose>
      </li>
     </ul>
     <ul class="none">
      <li>
       <h3>Keywords</h3>
       <xsl:choose>
        <xsl:when test="count(SECTION_PARAMETERS) > 0 or count(DEFAULT_KEYWORD) > 0 or count(KEYWORD) > 0">
         <xsl:if test="count(SECTION_PARAMETERS) > 0">
          <xsl:call-template name="list_keywords">
           <xsl:with-param name="element" select="SECTION_PARAMETERS"/>
          </xsl:call-template>
         </xsl:if>
         <xsl:if test="count(DEFAULT_KEYWORD) > 0">
          <xsl:call-template name="list_keywords">
           <xsl:with-param name="element" select="DEFAULT_KEYWORD"/>
          </xsl:call-template>
         </xsl:if>
         <xsl:if test="count(KEYWORD) > 0">
          <xsl:call-template name="list_keywords">
           <xsl:with-param name="element" select="KEYWORD"/>
          </xsl:call-template>
         </xsl:if>
        </xsl:when>
        <xsl:otherwise>
         <ul class="disc">
          <li>
           none
          </li>
         </ul>
        </xsl:otherwise>
       </xsl:choose>
      </li>
     </ul>
     <xsl:if test="count(SECTION_PARAMETERS) > 0 or count(DEFAULT_KEYWORD) > 0 or count(KEYWORD) > 0">
      <ul class="none">
       <li>
        <h3>Keyword descriptions</h3>
        <xsl:if test="count(SECTION_PARAMETERS) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="SECTION_PARAMETERS"/>
          <xsl:with-param name="root" select="$root"/>
          <xsl:with-param name="section_filename" select="$section_filename"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(DEFAULT_KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="DEFAULT_KEYWORD"/>
          <xsl:with-param name="root" select="$root"/>
          <xsl:with-param name="section_filename" select="$section_filename"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="KEYWORD"/>
          <xsl:with-param name="root" select="$root"/>
          <xsl:with-param name="section_filename" select="$section_filename"/>
         </xsl:call-template>
        </xsl:if>
       </li>
      </ul>
     </xsl:if>
     <xsl:call-template name="footer">
      <xsl:with-param name="root" select="$root"/>
     </xsl:call-template>
    </body>
   </html>
   </xsl:result-document>
   <xsl:call-template name="describe_sections">
    <xsl:with-param name="path" select="$local_path"/>
    <xsl:with-param name="root" select="concat($root,'../')"/>
   </xsl:call-template>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="list_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <ul class="disc">
  <xsl:for-each select="$element">
   <xsl:sort select="NAME[@type='default']"/>
   <xsl:if test="not(starts-with(NAME[@type='default'],'__'))">
    <li>
     <a href="#desc_{string(NAME[@type='default'])}" id="list_{string(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>
    </li>
   </xsl:if>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <xsl:param name="root"/>
 <xsl:param name="section_filename"/>
 <xsl:for-each select="$element">
  <xsl:sort select="NAME[@type='default']"/>
  <xsl:if test="not(starts-with(NAME[@type='default'],'__'))">
   <table class="default">
    <tr>
     <td class="l">
      <ul class="disc">
       <li>
        <a href="#list_{string(NAME[@type='default'])}" id="desc_{string(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>
       </li>
      </ul>
     </td>
     <td class="r">
     </td>
    </tr>
    <tr>
     <td class="l">
     </td>
     <td class="r">
      <xsl:if test="string-length(USAGE) > 0">
       <big class="uctt">
        <xsl:value-of select="USAGE"/>
       </big>
      </xsl:if>
     </td>
    </tr>
    <tr>
     <td class="l">
     </td>
     <td class="r">
      <em>
       <xsl:value-of select="DESCRIPTION"/>
      </em>
     </td>
    </tr>
    <tr>
     <td class="l">
     </td>
     <td class="r">
      This
      <xsl:if test="@required = 'yes'">
       required
      </xsl:if>
      <xsl:if test="@required = 'no'">
       optional
      </xsl:if>
      keyword can<xsl:if test="@repeats = 'no'">not</xsl:if> be repeated
      and it expects
      <xsl:if test="DATA_TYPE/N_VAR = -1">
       a list of <xsl:value-of select="DATA_TYPE/@kind"/>s.
      </xsl:if>
      <xsl:if test="DATA_TYPE/N_VAR = 1">
       precisely one <xsl:value-of select="DATA_TYPE/@kind"/>.
      </xsl:if>
      <xsl:if test="DATA_TYPE/N_VAR > 1">
       precisely <xsl:value-of select="DATA_TYPE/N_VAR"/>&#160;<xsl:value-of select="DATA_TYPE/@kind"/>s.
      </xsl:if>
     </td>
    </tr>
    <xsl:if test="LONE_KEYWORD_VALUE">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       The lone keyword behaves as a switch to
       <big class="uctt">
        <xsl:if test="LONE_KEYWORD_VALUE = 'T'">.TRUE.</xsl:if>
        <xsl:if test="LONE_KEYWORD_VALUE = 'F'">.FALSE.</xsl:if>
       </big>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="string-length(DEFAULT_VALUE) > 0">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       Default value<xsl:if test="DATA_TYPE/N_VAR > 1">s</xsl:if>:
       <big class="uctt">
        <xsl:choose>
         <xsl:when test="DATA_TYPE/@kind = 'logical'">
          <xsl:if test="DEFAULT_VALUE = 'T'">.TRUE.</xsl:if>
          <xsl:if test="DEFAULT_VALUE = 'F'">.FALSE.</xsl:if>
         </xsl:when>
         <xsl:otherwise>
          <xsl:value-of select="DEFAULT_VALUE"/>
         </xsl:otherwise>
        </xsl:choose>
       </big>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="string-length(DEFAULT_UNIT) > 0">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       Default unit:
       <big class="tt">
        [<xsl:value-of select="DEFAULT_UNIT"/>]
       </big>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="DATA_TYPE/ENUMERATION">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       List of valid keywords:
       <ul class="square">
        <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
         <xsl:sort select="NAME"/>
         <li>
          <dl>
           <dt>
            <big class="uctt">
             <xsl:value-of select="NAME"/>
            </big>
           </dt>
           <xsl:if test="string-length(DESCRIPTION) > 0">
            <dd>
             <em>
              <xsl:value-of select="DESCRIPTION"/>
             </em>
            </dd>
           </xsl:if>
          </dl>
         </li>
        </xsl:for-each>
       </ul>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="NAME[@type='alias']">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       Alias names for this keyword:
       <xsl:for-each select="NAME[@type='alias']">
        <xsl:sort select="NAME[@type='alias']"/>
        <xsl:choose>
         <xsl:when test="position() = last()">
          <xsl:value-of select="."/>
         </xsl:when>
         <xsl:otherwise>
          <xsl:value-of select="."/>,
         </xsl:otherwise>
        </xsl:choose>
       </xsl:for-each>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="count(REFERENCE) > 0">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       This keyword cites the following reference<xsl:if test="count(REFERENCE) > 1">s</xsl:if>:
       <xsl:for-each select="REFERENCE">
        <xsl:sort select="NAME"/>
        [<a href="{$root}references.html#reference_{string(NUMBER)}"><xsl:value-of select="NAME"/></a>]
       </xsl:for-each>
      </td>
     </tr>
    </xsl:if>
   </table>
  </xsl:if>
 </xsl:for-each>
</xsl:template>

<xsl:template name="link_section_path">
 <xsl:param name="string"/>
 <xsl:param name="separator"/>
 <xsl:param name="root"/>
 <xsl:variable name="string_before" select="substring-before($string,$separator)"/>
 <xsl:variable name="string_after" select="substring-after($string,$separator)"/>
 <a href="{concat($root,$string_before,'.html')}"><xsl:value-of select="$string_before"/></a> /
 <xsl:choose>
  <xsl:when test="contains($string_after,$separator)">
   <xsl:call-template name="link_section_path">
    <xsl:with-param name="string" select="$string_after"/>
    <xsl:with-param name="separator" select="$separator"/>
    <xsl:with-param name="root" select="substring-after($root,'../')"/>
   </xsl:call-template>
  </xsl:when>
  <xsl:otherwise>
   <a href="{concat($string_after,'.html')}"><xsl:value-of select="$string_after"/></a>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

</xsl:stylesheet>
