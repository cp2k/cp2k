<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="2.0">

<xsl:output method="html" indent="yes" name="html"/>

<xsl:template match="/CP2K_INPUT">
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
   <title>CP2K input reference</title>
   <style type="text/css">
    body {background-color: #ffffff}
    big.tt {font-family: monospace}
    big.uctt {font-family: monospace; text-transform: uppercase}
    ul.none {list-style-type: none}
    ul.disc {list-style-type: disc}
    ul.circle {list-style-type: circle}
   </style>
   <script src="collapsibleList.js" type="text/javascript" language="javascript1.2"/>
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
   <xsl:call-template name="describe_sections">
    <xsl:with-param name="path" select="'CP2K_INPUT'"/>
    <xsl:with-param name="root" select="'../'"/>
   </xsl:call-template>
   <hr/>
   <p>
    <a href="javascript:history.back()"><img src="l_arrow.gif"/></a>
    <a href="#top"><img src="u_arrow.gif"/></a>
   </p>
   <p>
    Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
   </p>
   <hr/>
  </body>
 </html>
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
    <a href="{$section_filename}"><xsl:value-of select="NAME"/></a>
   </li>
   <xsl:result-document href="{$section_filename}" method="html" indent="yes" format="html">
     <h3><a href="javascript:history.back()">Section &amp;<xsl:value-of select="NAME"/></a></h3>
    <xsl:if test="string-length(DESCRIPTION) > 0">
     <ul class="none">
      <li>
       <em>
        <xsl:value-of select="DESCRIPTION"/>
       </em>
      </li>
     </ul>
    </xsl:if>
    <xsl:if test="count(REFERENCE) > 0">
     <ul class="none">
      <li>
       References:
       <xsl:for-each select="REFERENCE">
        <xsl:sort select="NAME"/>
        [<a href="{$root}references.html#reference_{string(NUMBER)}"><xsl:value-of select="NAME"/></a>]
       </xsl:for-each>
      </li>
     </ul>
    </xsl:if>
    <ul class="none">
     <li>
      Section path:
      <big class="tt">
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
      This section can be repeated:
      <big class="uctt">
       <xsl:value-of select="@repeats"/>
      </big>
     </li>
    </ul>
    <ul class="none">
     <li>
      This section is required:
      <big class="uctt">
       <xsl:value-of select="@required"/>
      </big>
     </li>
    </ul>
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
        <ul class="none">
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
      <h3>Keywords:</h3>
      <xsl:choose>
       <xsl:when test="count(SECTION_PARAMETERS) > 0 or count(DEFAULT_KEYWORD) > 0 or count(KEYWORD) > 0">
        <xsl:if test="count(SECTION_PARAMETERS) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="SECTION_PARAMETERS"/>
          <xsl:with-param name="root" select="$root"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(DEFAULT_KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="element" select="DEFAULT_KEYWORD"/>
          <xsl:with-param name="root" select="$root"/>
         </xsl:call-template>
        </xsl:if>
        <xsl:if test="count(KEYWORD) > 0">
         <xsl:call-template name="describe_keywords">
          <xsl:with-param name="root" select="$root"/>
         </xsl:call-template>
        </xsl:if>
       </xsl:when>
       <xsl:otherwise>
        <ul class="none">
         <li>
          none
         </li>
        </ul>
       </xsl:otherwise>
      </xsl:choose>
     </li>
    </ul>
    <hr/>
    <p>
      <a href="javascript:history.back()"><img src="{$root}l_arrow.gif"/></a>
     <a href="#top"><img src="{$root}u_arrow.gif"/></a>
    </p>
    <p>
     Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
    </p>
    <hr/>
   </xsl:result-document>
   <xsl:call-template name="describe_sections">
    <xsl:with-param name="path" select="$local_path"/>
    <xsl:with-param name="root" select="concat($root,'../')"/>
   </xsl:call-template>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <xsl:param name="root"/>
 <ul class="none">
  <xsl:for-each select="$element">
   <xsl:sort select="NAME[@type='default']"/>
   <xsl:if test="not(starts-with(NAME[@type='default'],'__'))">
    <li>
     <dl>
      <dt>
       <u><xsl:value-of select="NAME[@type='default']"/></u>
       <xsl:if test="NAME[@type='alias']">
        (alias:
        <xsl:for-each select="NAME[@type='alias']">
         <xsl:sort select="NAME[@type='alias']"/>
         <xsl:choose>
          <xsl:when test="position() = last()">
           <xsl:value-of select="."/>)
          </xsl:when>
          <xsl:otherwise>
           <xsl:value-of select="."/>,
          </xsl:otherwise>
         </xsl:choose>
        </xsl:for-each>
       </xsl:if>
      </dt>
      <dd>
       <p>
        <em>
         <xsl:value-of select="DESCRIPTION"/>
        </em>
       </p>
      </dd>
      <xsl:if test="count(REFERENCE) > 0">
       <dd>
        <p>
         References:
         <xsl:for-each select="REFERENCE">
          <xsl:sort select="NAME"/>
          [<a href="{$root}references.html#reference_{string(NUMBER)}"><xsl:value-of select="NAME"/></a>]
         </xsl:for-each>
        </p>
       </dd>
      </xsl:if>
      <dd>
       <p>
        This keyword can be repeated:
        <big class="uctt">
         <xsl:value-of select="@repeats"/>
        </big>
       </p>
      </dd>
      <dd>
       <p>
        This keyword is required:
        <big class="uctt">
         <xsl:value-of select="@required"/>
        </big>
       </p>
      </dd>
      <dd>
       <p>
        Data type: <big class="uctt"><xsl:value-of select="DATA_TYPE/@kind"/></big>
       </p>
       <xsl:if test="DATA_TYPE/ENUMERATION">
        <p>
         List of valid keys:
        </p>
        <ul class="none">
         <li>
          <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
           <xsl:sort select="NAME"/>
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
          </xsl:for-each>
         </li>
        </ul>
       </xsl:if>
      </dd>
      <xsl:if test="string-length(USAGE) > 0">
       <dd>
        <p>
         Usage example:
         <big class="uctt">
          <xsl:value-of select="USAGE"/>
         </big>
        </p>
       </dd>
      </xsl:if>
      <xsl:if test="string-length(DEFAULT_VALUE) > 0">
       <dd>
        <p>
         Default value:
         <big class="uctt">
          <xsl:value-of select="DEFAULT_VALUE"/>
         </big>
        </p>
       </dd>
      </xsl:if>
      <xsl:if test="string-length(DEFAULT_UNIT) > 0">
       <dd>
        <p>
         Default unit:
         <big class="uctt">
          <xsl:value-of select="DEFAULT_UNIT"/>
         </big>
        </p>
       </dd>
      </xsl:if>
     </dl>
    </li>
   </xsl:if>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="link_section_path">
 <xsl:param name="string"/>
 <xsl:param name="separator"/>
 <xsl:param name="root"/>
 <xsl:variable name="string_before" select="substring-before($string,$separator)"/>
 <xsl:variable name="string_after" select="substring-after($string,$separator)"/>
 <a href="{concat($root,$string_before,'.html')}"><xsl:value-of select="concat($string_before,'/')"/></a>
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
