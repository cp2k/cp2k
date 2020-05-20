<?xml version="1.0" encoding="UTF-8"?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">

<xsl:output doctype-public="html" doctype-system="html" indent="yes" method="html" name="html"/>

<xsl:param name="add_edit_links" select="'no'"/>
<xsl:param name="version" select="'trunk'"/>
<xsl:param name="release_path" select="'/trunk'"/>

<xsl:template match="/CP2K_INPUT">
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
   <title>CP2K input reference</title>
   <xsl:call-template name="head">
    <xsl:with-param name="root" select="''"/>
   </xsl:call-template>
  </head>
  <body>
   <xsl:call-template name="write_generate_manual_howto"/>
   <xsl:call-template name="write_html_tables"/>
   <h1 align="center">CP2K input reference</h1>
   <h2>Version information</h2>
   <p>
    This HTML manual refers to
    <a href="https://github.com/cp2k/cp2k/tree/{substring-after(/CP2K_INPUT/COMPILE_REVISION,':')}" target="_blank">
     <xsl:value-of select="/CP2K_INPUT/CP2K_VERSION"/> (<xsl:value-of select="/CP2K_INPUT/COMPILE_REVISION"/>)
    </a>
    and was generated automatically from a CP2K executable compiled on <xsl:value-of select="COMPILE_DATE"/> using the
    <code>--xml</code> command line option (see <a href="generate_manual_howto.html">how to generate this manual</a>).
    Thus the manual describes exactly this version of the code.
   </p>
   <p>
    <xsl:call-template name="searchform"/>
   </p>
   <h2>Journal papers</h2>
   <p>
    <a href="references.html">List of references</a> cited in the CP2K input manual.
   </p>
   <h2>CP2K units</h2>
   <p>
    <a href="units.html">Available units of measurement</a> which can be used in the CP2K input for keyword values.
    The units can be used to specify an alternative unit of measurement for keyword values, for which a
    default unit has been explicitly defined.
    They are specified in the form of a string <code>[unit]</code> as a separate token before any value.
    A unit specification is applied to all values on the same line following that specification, unless explicitly
    overwritten by a another unit specification.
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
     recalled with a <i>${VAR}</i> (or  <i>$VAR</i>) statement. There can be only one @SET statement per line.
     A variable name can only contain alphanumerical characters plus the underscore and must not start with a number.
    </dd>
    <dt><b>${VAR}</b> or <b>$VAR</b></dt>
    <dd>
     Expand the variable <i>VAR</i>. The text <i>${VAR}</i> (or  <i>$VAR</i>) is replaced
     with the value assigned to <i>VAR</i> in the last @SET directive.
     There can be multiple variable statements per line. The expansion process is repeated
     until no more variables are found.
    </dd>
    <dt><b>${VAR-DEFAULT}</b></dt>
    <dd>
     Expand the variable <i>VAR</i> if defined. The text <i>${VAR-DEFAULT}</i> is replaced
     with the value assigned to <i>VAR</i> in the last @SET directive if defined,
     otherwise with the value specified after the dash <i>-</i> (here: <i>DEFAULT</i>).
    </dd>
    <dt><b>@IF / @ENDIF</b></dt>
    <dd>
     Conditional block. The text from the @IF line up to the next line with a valid
     @ENDIF is skipped, if the expression following @IF resolves to <i>false</i>.
     Available expressions are lexical comparisons for equality '==' or inequality '/='.
     If none of the two operators are found, a '0' or whitespace resolves to <i>false</i>
     while any text resolves to <i>true</i>. @IF/@ENDIF blocks cannot be nested and
     cannot span across files. There can be only one test (== or /=) per @IF statement.
    </dd>
   </dl>
   <h2>Input structure</h2>
    All sections that can be part of a CP2K input file are shown here with their allowed nestings.
    A detailed description of each section and its keywords can be obtained by clicking on the
    section links. The links in the detailed descriptions switch back to the corresponding index
    entries. In this way a toggling between the index and the detailed description is feasible.
   <h2>Index of all input sections</h2>
   <ul class="noscript"> 
    <li>
     Double click on [&#8722;] / [+] to shrink or expand all subtrees contained in a section
    </li>
    <li>
     Single click on [&#8722;] / [+] to shrink or expand the top level subtree contained in a section
    </li>
   </ul>
   <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
   <ul class="index">
    <li>
     <a href="CP2K_INPUT.html" id="CP2K_INPUT.html">CP2K_INPUT</a>
    </li>
    <xsl:call-template name="describe_sections">
     <xsl:with-param name="path" select="'CP2K_INPUT'"/>
     <xsl:with-param name="root" select="'../'"/>
    </xsl:call-template>
   </ul>
   <xsl:call-template name="footer">
    <xsl:with-param name="root" select="''"/>
   </xsl:call-template>
  </body>
  <xsl:result-document href="CP2K_INPUT.html" method="html" indent="yes" format="html">
   <html>
    <head>
     <title>CP2K input file</title>
     <xsl:call-template name="head">
      <xsl:with-param name="root" select="''"/>
     </xsl:call-template>
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
       <code class="upper">
        <a href="CP2K_INPUT.html">CP2K_INPUT</a>
       </code>
      </li>
     </ul>
     <ul class="none">
      <li>
       This section cannot be repeated.
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
 <xsl:param name="root"/>
 <xsl:variable name="description">
  <xsl:choose>
   <xsl:when test="string-length(DESCRIPTION) > 0">
    <xsl:value-of disable-output-escaping="yes" select="string(DESCRIPTION)"/>
   </xsl:when>
   <xsl:otherwise>
    <xsl:value-of disable-output-escaping="yes" select="'CP2K input reference'"/>
   </xsl:otherwise>
  </xsl:choose>
 </xsl:variable>
 <meta name="language" content="en"/>
 <meta name="copyright" content="2000 - {year-from-dateTime(current-dateTime())} CP2K developers group"/>
 <meta name="description" content="{$description}"/>
 <meta name="keywords" content="scientific,computing,chemistry,physics,documentation,help,manual,Fortran,parallel,molecular dynamics,MD,density functional theory,DFT,electronic structure,linear scaling,force field,Quickstep,GPW,GAPW,FIST,QM,MM"/>
 <link rel="shortcut icon" href="{$root}favicon.png" type="image/png"/>
 <style type="text/css">
  a {text-decoration: none;}
  body {background-color: #ffffff;}
  code.upper {text-transform: uppercase;}
  li {margin-left: 0em; padding-left: 0em; text-indent: 0em;}
  p.uctt {font-family: monospace; text-transform: uppercase;}
  table.default {table-layout: fixed; width: 100%;}
  td.l {width: 25%;}
  td.r {width: 75%;}
  ul.circle {list-style-type: circle;}
  ul.disc {list-style-type: disc;}
  ul.index {list-style-type: none; margin-left: 0em; padding-left: 1.8em; text-indent: 0em;}
  ul.none {list-style-type: none;}
  ul.noscript {list-style-type: disc;}
  ul.square {list-style-type: square;}
  .button {font-family: monospace; font-size: 100%; cursor: pointer;}
  #html_table
  {
   border: 1px solid #000000;
   border-collapse: collapse;
   margin-left: 25px;
   padding: 6px;
   text-align: left;
   vertical-align: middle;
  }
 </style>
 <noscript>
  <style>
   ul.index {list-style-type: disc; margin-left: 0px; padding-left: 1.8em; text-indent: 0px}
   ul.noscript {display: none}
  </style>
 </noscript>
 <script language="javascript" type="text/javascript" src="{$root}toggle_folding.js"></script>
</xsl:template>

<xsl:template name="header">
 <xsl:param name="root"/>
 <table class="default" summary="header">
  <tr>
   <td align="left">
    Back to the <a href="{$root}index.html">main page</a> of this manual
   </td>
   <td align="center">
    Input reference of <xsl:value-of select="/CP2K_INPUT/CP2K_VERSION"/>
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
 <table class="default" summary="footer">
  <tr>
   <td align="left">
    Back to the <a href="{$root}index.html">main page</a> of this manual or the <a href="https://www.cp2k.org">CP2K home page</a>
   </td>
   <td align="right">
    (Last update:
    <xsl:value-of select="day-from-dateTime(current-dateTime())"/>.<xsl:value-of select="month-from-dateTime(current-dateTime())"/>.<xsl:value-of select="year-from-dateTime(current-dateTime())"/>)
   </td>
  </tr>
 </table>
</xsl:template>

<xsl:template name="searchform">
 <form method="get" action="https://www.google.com/search">
  <input type="text" name="q" maxlength="255"/>
  <input type="submit" value="Search this manual (Google)"/>
  <input type="hidden" name="domains" value="https://manual.cp2k.org/{$version}/"/>
  <input type="radio" style="visibility:hidden" name="sitesearch" value="https://manual.cp2k.org/{$version}/" checked="checked"/>
 </form>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:param name="path"/>
 <xsl:param name="root"/>
 <ul class="index">
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
     <title><xsl:value-of select="replace($local_path,'/',' / ')" /></title>
     <xsl:call-template name="head">
      <xsl:with-param name="root" select="$root"/>
     </xsl:call-template>
    </head>
    <body>
     <xsl:call-template name="header">
      <xsl:with-param name="root" select="$root"/>
     </xsl:call-template>
     <h2><a href="{$root}index.html#{$section_filename}">Section <xsl:value-of select="NAME"/></a></h2>
     <ul class="none">
      <li>
       <em>
        <xsl:if test="string-length(DESCRIPTION) > 0">
         <xsl:value-of disable-output-escaping="yes" select="DESCRIPTION"/>
        </xsl:if>
        <xsl:if test="string-length(DESCRIPTION) = 0">
         Without description, yet.
        </xsl:if>
       </em>
       <xsl:call-template name="link_edit_text">
        <xsl:with-param name="location" select="string(LOCATION)"/>
       </xsl:call-template>
      </li>
     </ul>
     <ul class="none">
      <li>
       Section path:
       <code class="upper">
        <xsl:call-template name="link_section_path">
         <xsl:with-param name="string" select="$local_path"/>
         <xsl:with-param name="separator" select="'/'"/>
         <xsl:with-param name="root" select="$root"/>
        </xsl:call-template>
       </code>
      </li>
     </ul>
     <ul class="none">
      <li>
       This section can<xsl:if test="@repeats = 'no'">not</xsl:if> be repeated.
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
          <xsl:with-param name="local_path" select="$local_path"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(DEFAULT_KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="DEFAULT_KEYWORD"/>
          <xsl:with-param name="root" select="$root"/>
          <xsl:with-param name="section_filename" select="$section_filename"/>
          <xsl:with-param name="local_path" select="$local_path"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="KEYWORD"/>
          <xsl:with-param name="root" select="$root"/>
          <xsl:with-param name="section_filename" select="$section_filename"/>
          <xsl:with-param name="local_path" select="$local_path"/>
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
   <xsl:if test="not(starts-with(NAME[@type='default'],'__CONTROL'))">
    <li>
     <a href="#{string(NAME[@type='default'])}" id="list_{string(NAME[@type='default'])}">
     <xsl:call-template name="keyword_name">
      <xsl:with-param name="element" select="$element"/>
     </xsl:call-template>
     </a>
    </li>
   </xsl:if>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="keyword_name">
 <xsl:param name="element" select="KEYWORD"/>
 <xsl:choose>
  <xsl:when test="@removed = 'yes' or string-length(DEPRECATION_NOTICE) > 0 ">
   <span style="text-decoration: line-through">
    <xsl:value-of select="NAME[@type='default']"/>
   </span>
  </xsl:when>
  <xsl:otherwise>
   <xsl:value-of select="NAME[@type='default']"/>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <xsl:param name="root"/>
 <xsl:param name="section_filename"/>
 <xsl:param name="local_path"/>
 <xsl:for-each select="$element">
  <xsl:sort select="NAME[@type='default']"/>
  <xsl:if test="not(starts-with(NAME[@type='default'],'__CONTROL'))">
  <xsl:variable name="keyword_path" select="concat($local_path,'/',string(NAME[@type='default']))"/>
   <table class="default" summary="keyword_description">
    <tr>
     <td class="l">
      <ul class="disc">
       <li>
        <a id="desc_{string(NAME[@type='default'])}"></a>
        <a href="#list_{string(NAME[@type='default'])}" id="{string(NAME[@type='default'])}">
         <xsl:call-template name="keyword_name">
          <xsl:with-param name="element" select="$element"/>
         </xsl:call-template>
        </a>
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
      <xsl:choose>
       <xsl:when test="NAME[@type='default'] = 'DEFAULT_KEYWORD'">
        <code><xsl:value-of disable-output-escaping="yes" select="USAGE"/></code>
       </xsl:when>
       <xsl:otherwise>
        <xsl:variable name="vartype" select="concat(upper-case(substring(DATA_TYPE/@kind,1,1)),substring(DATA_TYPE/@kind,2))"/>
        <xsl:if test="DATA_TYPE/N_VAR > 0">
         <xsl:choose>
          <xsl:when test="NAME[@type='default'] = 'SECTION_PARAMETERS'">
           <code class="upper">&amp;<xsl:value-of disable-output-escaping="yes" select="../NAME"/></code>
          </xsl:when>
          <xsl:otherwise>
           <code class="upper"><xsl:value-of disable-output-escaping="yes" select="NAME[@type='default']"/></code>
          </xsl:otherwise>
         </xsl:choose>
         <xsl:call-template name="repeat">
          <xsl:with-param name="ivar" select="1"/>
          <xsl:with-param name="nvar" select="DATA_TYPE/N_VAR"/>
          <xsl:with-param name="vartype" select="$vartype"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="DATA_TYPE/N_VAR = -1">
         <code class="upper"><xsl:value-of disable-output-escaping="yes" select="NAME[@type='default']"/></code>
         <code>&#160;{<xsl:value-of select="$vartype"/>}&#160;...</code>
         <xsl:if test="contains(upper-case(NAME[@type='default']),'LIST')">
          <code>&#160;or&#160;a&#160;range&#160;{<xsl:value-of select="$vartype"/>}..{<xsl:value-of select="$vartype"/>}</code>
         </xsl:if>
        </xsl:if>
       </xsl:otherwise>
      </xsl:choose>
     </td>
    </tr>
    <tr>
     <td class="l">
     </td>
     <td class="r">
      <em>
       <xsl:if test="@removed = 'yes'">
        <p style="font-weight:bold">
         This keyword was removed.
         <xsl:value-of select="DEPRECATION_NOTICE"/>
        </p>
       </xsl:if>
       <xsl:if test="@removed = 'no' and string-length(DEPRECATION_NOTICE) > 0">
        <p style="font-weight:bold">
         This keyword is deprecated.
         <xsl:value-of select="DEPRECATION_NOTICE"/>
        </p>
       </xsl:if>
       <xsl:if test="string-length(DESCRIPTION) > 0">
        <xsl:value-of disable-output-escaping="yes" select="DESCRIPTION"/>
       </xsl:if>
       <xsl:if test="string-length(DESCRIPTION) = 0">
        Without description, yet.
       </xsl:if>
      </em>
      <xsl:call-template name="link_edit_text">
       <xsl:with-param name="location" select="string(LOCATION)"/>
      </xsl:call-template>
     </td>
    </tr>
    <tr>
     <td class="l">
     </td>
     <td class="r">
      This keyword can<xsl:if test="@repeats = 'no'">not</xsl:if> be repeated
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
       <xsl:if test="DATA_TYPE/@kind = 'logical'">
        The lone keyword behaves as a switch to
        <code class="upper">
         <xsl:if test="LONE_KEYWORD_VALUE = 'T'">.TRUE.</xsl:if>
         <xsl:if test="LONE_KEYWORD_VALUE = 'F'">.FALSE.</xsl:if>
        </code>
       </xsl:if>
       <xsl:if test="DATA_TYPE/@kind = 'integer' or DATA_TYPE/@kind = 'real'">
        The lone keyword defaults to
        <code class="upper">
         <xsl:value-of select="LONE_KEYWORD_VALUE"/>
        </code>
       </xsl:if>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="string-length(DEFAULT_VALUE) > 0">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       Default value<xsl:if test="DATA_TYPE/N_VAR > 1">s</xsl:if>:
       <code class="upper">
        <xsl:choose>
         <xsl:when test="DATA_TYPE/@kind = 'logical'">
          <xsl:if test="DEFAULT_VALUE = 'T'">.TRUE.</xsl:if>
          <xsl:if test="DEFAULT_VALUE = 'F'">.FALSE.</xsl:if>
         </xsl:when>
         <xsl:otherwise>
          <xsl:value-of select="DEFAULT_VALUE"/>
         </xsl:otherwise>
        </xsl:choose>
       </code>
      </td>
     </tr>
    </xsl:if>
    <xsl:if test="string-length(DEFAULT_UNIT) > 0">
     <tr>
      <td class="l">
      </td>
      <td class="r">
       <a href="{$root}units.html">Default unit:</a>
       <code>
        [<xsl:value-of select="DEFAULT_UNIT"/>]
       </code>
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
            <code class="upper">
             <xsl:value-of select="NAME"/>
            </code>
           </dt>
           <dd>
            <em>
             <xsl:if test="string-length(DESCRIPTION) > 0">
              <xsl:value-of disable-output-escaping="yes" select="DESCRIPTION"/>
             </xsl:if>
             <xsl:if test="string-length(DESCRIPTION) = 0">
              Without description, yet.
             </xsl:if>
            </em>
           </dd>
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

<xsl:template name="write_generate_manual_howto">
 <xsl:result-document href="generate_manual_howto.html" method="html" indent="yes" format="html">
  <html>
   <head>
    <title>How to generate the CP2K input reference manual</title>
    <xsl:call-template name="head">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
   </head>
   <body>
    <xsl:call-template name="header">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
    <h2 align="center">How to generate the CP2K input reference manual</h2>
    <ul class="disc">
     <li>
      Prepare a CP2K executable. It does not matter which type of CP2K executable (e.g. psmp) you are using,
      since running a CP2K executable with the flag <code>--xml</code> like:
      <ul class="none">
       <li><code>cp2k.psmp --xml</code></li>
      </ul>
      will generate a file named <code>cp2k_input.xml</code> with an XML dump of the entire CP2K input structure
      in the same directory describing all features of the actual CP2K executable.<br/>
      That XML output in <code>cp2k_input.xml</code> has to be transformed to HTML using an XML 2.0 compliant
      translator like SAXON. Preferentially, you may want to employ the platform independent Java version of SAXON, which
      can be downloaded from <a href="https://sourceforge.net/projects/saxon">https://sourceforge.net/projects/saxon</a>.<br/>
     </li>
     <li>
      Alternatively, you can also run the script
      <ul class="none">
       <li><code>update_manual local psmp</code></li>
      </ul>
      in the folder cp2k/tools/manual/ which will perform all steps using the CP2K executable cp2k/exe/local/cp2k.psmp.
     </li>
     <li>
      Launch your favorite web browser and load the generated <code>index.html</code> file, e.g. in cp2k/tools/manual/local/psmp
      from the previous step.
     </li>
     <li>
      Happy browsing!
     </li>
    </ul>
    <xsl:call-template name="footer">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
   </body>
  </html>
 </xsl:result-document>
</xsl:template>

<xsl:template name="link_edit_text">
 <xsl:param name="location"/>
 <xsl:if test="$add_edit_links = 'yes'">
 <xsl:variable name="link_path" select="concat(substring-before($location, ':'), '#L', substring-after($location, ':'))"/>
 <xsl:variable name="link_url" select="concat('https://github.com/cp2k/cp2k/blob/master/src/', $link_path)"/>
  <span style="font-size: small;">
   &#160;[<a title="See corresponding source code location." href="{$link_url}">Edit on GitHub</a>]
  </span>
 </xsl:if>
</xsl:template>

<xsl:template name="write_html_tables">
 <xsl:result-document href="html_tables.html" method="html" indent="yes" format="html">
  <html>
   <head>
    <title>Supported HTML entities and tags</title>
    <xsl:call-template name="head">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
   </head>
   <body>
    <xsl:call-template name="header">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
    <h2 align="center">List of HTML tags supported by CP2K</h2>
    <p>The following HTML tags are supported by CP2K within the input descriptions:</p>
    <ul class="disc">
     <xsl:for-each select="CP2K_HTML/TAG">
      <li>
       <xsl:value-of disable-output-escaping="no" select="NAME"/>
      </li>
     </xsl:for-each>
    </ul>
    <h2 align="center">Table of HTML entities supported by CP2K</h2>
    <p>The following HTML entities are supported by CP2K within the input descriptions:</p>
    <table id="html_table">
     <tr id="html_table">
      <th id="html_table">Name</th>
      <th id="html_table">Code</th>
      <th id="html_table">Character</th>
     </tr>
     <xsl:for-each select="CP2K_HTML/ENTITY">
      <tr id="html_table">
       <td id="html_table"><xsl:value-of disable-output-escaping="no" select="NAME"/></td>
       <td id="html_table"><xsl:value-of disable-output-escaping="no" select="CODE"/></td>
       <td id="html_table"><xsl:value-of disable-output-escaping="yes" select="CODE"/></td>
      </tr>
     </xsl:for-each>
    </table>
    <xsl:call-template name="footer">
     <xsl:with-param name="root" select="''"/>
    </xsl:call-template>
   </body>
  </html>
 </xsl:result-document>
</xsl:template>

<xsl:template name="repeat">
 <xsl:param name="ivar"/>
 <xsl:param name="nvar"/>
 <xsl:param name="vartype"/>
 <code>&#160;{<xsl:value-of select="concat(upper-case(substring($vartype,1,1)),substring($vartype,2))"/>}</code>
 <xsl:if test="not($ivar = $nvar)">
  <xsl:call-template name="repeat">
   <xsl:with-param name="ivar" select="$ivar + 1"/>
   <xsl:with-param name="nvar" select="$nvar"/>
   <xsl:with-param name="vartype" select="$vartype"/>
  </xsl:call-template>
 </xsl:if>
</xsl:template>

</xsl:stylesheet>
